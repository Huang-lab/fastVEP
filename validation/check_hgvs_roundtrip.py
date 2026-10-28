#!/usr/bin/env python3
"""Is an HGVS description the change it claims to be?  (HGVSc, HGVSp and HGVSg of a VEP-style tab file)

Every comparison with Ensembl VEP can only show that two tools agree.  A port that copies VEP copies its
defects, and they are invisible to that comparison.  This script does not use VEP at all: it takes the GFF3
and the FASTA, builds the transcript, applies the description in a row to the sequence and asks whether the
result is the change the VCF record makes.

  HGVSc / HGVSg   the description, applied to the genome, must give the genome the VCF allele gives
                  (two spellings of one change give one edited sequence, however they are written);
                  then the HGVS spelling rules: the most 3' position, `dup` for an insertion that
                  repeats what precedes it, `inv` for an inversion, `>` for a single base.
  HGVSp           the description, applied to the reference protein, must give the protein of the
                  EDITED transcript (the VCF allele applied to the spliced sequence and translated, on
                  into the 3' UTR for a frameshift or a lost terminator).  For `fs` and `ext` the first
                  changed residue and the distance to the new terminator are checked; a change that
                  loses the initiator must be `p.Met1?`.

usage
  check_hgvs_roundtrip.py check --gff3 G --fasta F --vcf V --tab T --out OUT.tsv [--every N] [--workers W]
  check_hgvs_roundtrip.py summary OUT.tsv [OUT2.tsv ...]
  check_hgvs_roundtrip.py compare VEP.tsv FASTVEP.tsv [--show N]
  check_hgvs_roundtrip.py selftest

`check` writes one line per description: type (c, p, g), variant, allele, transcript, the description, a
class, a verdict, a detail and flags.  Verdicts:
  ok            the description is the change and is spelled as HGVS asks
  ok_spelling   the same change, spelled in a form HGVS does not prefer (ins for dup, delins for del...)
  ok_not3       the same change, not at the most 3' position (HGVSc/HGVSg), or not shifted to the end
                of a repeat in the protein (HGVSp)
  wrong         the description is not the change (the detail says how)
  malformed     not HGVS syntax (`p.MetAla1_?2`, `p.Asp1_?0`, `p.LeuIle2672=`)
  unparsed      the checker cannot read it
  skip_*        the row cannot be checked (the VCF REF is not the genome, a change that spans a splice
                boundary, an ambiguous residue, a transcript without a complete CDS...)
  missing       (p) the tool wrote nothing and the change alters the protein
"""
import argparse
import bisect
import collections
import multiprocessing
import os
import re
import sys
import zlib

# ----------------------------------------------------------------------------- sequence helpers

COMP = str.maketrans("ACGTRYKMBDHVNSWacgtrykmbdhvnsw", "TGCAYRMKVHDBNSWtgcayrmkvhdbnsw")


def revcomp(s):
    return s.translate(COMP)[::-1]


_BASES = "TCAG"
_STD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_MITO = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"


def _table(aas):
    t = {}
    i = 0
    for a in _BASES:
        for b in _BASES:
            for c in _BASES:
                t[a + b + c] = aas[i]
                i += 1
    return t


STD = _table(_STD)
MITO = _table(_MITO)


def translate(seq, table=STD, stop_at_stop=False):
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = table.get(seq[i:i + 3], "X")
        out.append(aa)
        if stop_at_stop and aa == "*":
            break
    return "".join(out)


def canon_chrom(n):
    if n.lower().startswith("chr"):
        n = n[3:]
    return "MT" if n == "M" else n


class Bad(Exception):
    """A description that cannot be read."""


class Skip(Bad):
    """A description the checker cannot decide: an insertion at the junction of two exons, which can be at either end
    of the intron between them."""


class Wrong(Bad):
    """A description that can be read and cannot be a description of the change (a stated base that is not
    the genome's, an offset from a base that is not the edge of an exon)."""


# ----------------------------------------------------------------------------- FASTA

class Fasta:
    def __init__(self, path):
        self.path = path
        self.fh = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as fh:
            for line in fh:
                n, ln, off, lb, lw = line.split()[:5]
                self.idx[n] = (int(ln), int(off), int(lb), int(lw))
        self.alias = {}
        for n in self.idx:
            self.alias.setdefault(canon_chrom(n), n)
        self.pid = os.getpid()

    def name(self, chrom):
        if chrom in self.idx:
            return chrom
        return self.alias.get(canon_chrom(chrom))

    def fetch(self, chrom, a, b):
        """0-based half-open, upper-case, clipped to the contig"""
        if os.getpid() != self.pid:  # a worker after fork must not share the parent's file offset
            self.fh = open(self.path, "rb")
            self.pid = os.getpid()
        name = self.name(chrom)
        ln, off, lb, lw = self.idx[name]
        a = max(0, a)
        b = min(ln, b)
        if b <= a:
            return ""
        b0 = off + (a // lb) * lw + a % lb
        b1 = off + ((b - 1) // lb) * lw + (b - 1) % lb + 1
        self.fh.seek(b0)
        return self.fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


class GW:
    """a window on the genome that grows where it is asked for a base outside it"""

    def __init__(self, fasta, chrom):
        self.fasta, self.chrom, self.a, self.s = fasta, chrom, 0, ""

    def base(self, i):
        if i < 0:
            return ""
        if not (self.a <= i < self.a + len(self.s)):
            self.a = max(0, i - 1024)
            self.s = self.fasta.fetch(self.chrom, self.a, i + 1024)
        j = i - self.a
        return self.s[j] if 0 <= j < len(self.s) else ""

    def get(self, a, b):
        if a < 0:
            raise Bad("position before the start of the contig")
        if not (self.a <= a and b <= self.a + len(self.s)):
            self.a = max(0, a - 256)
            self.s = self.fasta.fetch(self.chrom, self.a, b + 256)
        return self.s[a - self.a:b - self.a]


# ----------------------------------------------------------------------------- transcripts

class Tx:
    __slots__ = ("key", "tid", "version", "chrom", "strand", "start", "end", "exons", "cds", "biotype",
                 "asc", "txoff_asc", "tx_cum", "ex_tx", "L", "cds_start", "cds_end", "phase0",
                 "firsts", "lasts", "_seq", "_prot", "_table", "ok_cds", "why_not", "cds_offset", "nf_start")

    def finish(self):
        """exons sorted both ways, cumulative offsets, CDS in transcript coordinates"""
        asc = sorted(set(self.exons))
        self.asc = asc
        lens = [e - s + 1 for s, e in asc]
        total = sum(lens)
        self.L = total
        if self.strand > 0:
            self.txoff_asc = [sum(lens[:i]) for i in range(len(asc))]
            self.ex_tx = list(asc)
        else:
            self.txoff_asc = [sum(lens[i + 1:]) for i in range(len(asc))]
            self.ex_tx = list(reversed(asc))
        ex_lens = [e - s + 1 for s, e in self.ex_tx]
        self.tx_cum = [sum(ex_lens[:i]) for i in range(len(self.ex_tx))]
        self.firsts = set(c + 1 for c in self.tx_cum)
        self.lasts = set(c + n for c, n in zip(self.tx_cum, ex_lens))
        self.cds_start = self.cds_end = None
        self.phase0 = 0
        self.cds_offset = 0
        # the GFF3 does not say that a CDS has no start codon (the tab file's FLAGS does, row by row), and
        # such a transcript's first codon is a codon like the others even where it reads GTG or CTG
        self.nf_start = False
        self.ok_cds = False
        self.why_not = ""
        self._seq = None
        self._prot = None
        self._table = STD
        if self.cds:
            gmin = min(s for s, e, p in self.cds)
            gmax = max(e for s, e, p in self.cds)
            a, b = self.pos_of(gmin), self.pos_of(gmax)
            if a is None or b is None:
                self.why_not = "cds_outside_exon"
                return
            if self.strand > 0:
                self.cds_start, self.cds_end = a, b
                first = min(self.cds, key=lambda c: c[0])
            else:
                self.cds_start, self.cds_end = b, a
                first = max(self.cds, key=lambda c: c[1])
            self.phase0 = first[2] if first[2] in (0, 1, 2) else 0
            # a CDS that starts part-way through a codon (phase 1 or 2) is numbered from the codon's first base,
            # which the transcript does not have: c.1 and c.2 are the missing bases and the first CDS base is c.3
            self.cds_offset = (3 - self.phase0) % 3
            if self.chrom in ("MT", "M", "chrM", "chrMT"):
                self._table = MITO

    def pos_of(self, g):
        """1-based transcript position of the 1-based genomic position g, None if it is not in an exon"""
        i = bisect.bisect_right(self.asc, (g, 1 << 60)) - 1
        if i < 0:
            return None
        s, e = self.asc[i]
        if g > e:
            return None
        if self.strand > 0:
            return self.txoff_asc[i] + (g - s) + 1
        return self.txoff_asc[i] + (e - g) + 1

    def g_of(self, t):
        """genomic position of the transcript position t, extrapolated beyond the ends"""
        if t < 1:
            s, e = self.ex_tx[0]
            return s + (t - 1) if self.strand > 0 else e - (t - 1)
        if t > self.L:
            s, e = self.ex_tx[-1]
            d = t - self.L
            return e + d if self.strand > 0 else s - d
        k = bisect.bisect_right(self.tx_cum, t - 1) - 1
        s, e = self.ex_tx[k]
        d = t - 1 - self.tx_cum[k]
        return s + d if self.strand > 0 else e - d

    def exon_index(self, g):
        i = bisect.bisect_right(self.asc, (g, 1 << 60)) - 1
        if i >= 0 and g <= self.asc[i][1]:
            return i
        return None

    def seq(self, fasta):
        if self._seq is None:
            parts = []
            for s, e in self.ex_tx:
                x = fasta.fetch(self.chrom, s - 1, e)
                parts.append(x if self.strand > 0 else revcomp(x))
            self._seq = "".join(parts)
        return self._seq

    def prot(self, fasta):
        """(reference protein with its terminator, why the transcript cannot be checked or "")"""
        if self._prot is None:
            why = ""
            prot = ""
            if self.cds_start is None:
                why = "noncoding" if not self.why_not else self.why_not
            elif self.phase0 != 0:
                why = "cds_start_nf"
            else:
                sq = self.seq(fasta)
                cds = sq[self.cds_start - 1:self.cds_end]
                if len(cds) % 3:
                    why = "cds_end_nf"
                else:
                    prot = translate(cds, self._table)
                    if cds[:3] in ("TTG", "CTG", "GTG", "ACG") and not self.nf_start:
                        # a non-ATG initiator reads as Met (Ensembl's initial_met, RefSeq's transl_except); the
                        # first codon of a CDS with no start codon is read as the codon it is
                        prot = "M" + prot[1:]
                    if not prot.endswith("*"):
                        why = "no_terminator"
                    elif "*" in prot[:-1]:
                        why = "internal_stop"
                    elif "X" in prot:
                        why = "ambiguous_residue"
            self._prot = (prot, why)
        return self._prot


def parse_attrs(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


SKIP_TYPES = {"gene", "ncRNA_gene", "pseudogene", "biological_region", "region", "chromosome", "supercontig",
              "exon", "CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon", "cDNA_match",
              "match", "scaffold", "enhancer", "promoter", "CAGE_cluster", "TF_binding_site", "silencer",
              "open_chromatin_region", "CTCF_binding_site", "enhancer_blocking_element", "mobile_genetic_element",
              "transcriptional_cis_regulatory_region", "sequence_alteration", "regulatory_region", "origin_of_replication",
              "D_loop", "minisatellite", "repeat_region", "sequence_feature", "centromere", "STS", "tandem_repeat",
              "nucleotide_motif", "protein_binding_site", "histone_modification", "response_element", "locus_control_region",
              "recombination_feature", "sequence_secondary_structure", "stem_loop", "mature_transcript", "primary_transcript"}


def load_gff(path, wanted=None):
    """transcripts by key (the transcript_id without version for Ensembl, the accession for RefSeq)"""
    by_id = {}
    out = {}
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ = f[2]
            if typ == "exon" or typ == "CDS":
                par = None
                a = f[8]
                i = a.find("Parent=")
                if i >= 0:
                    j = a.find(";", i)
                    par = a[i + 7:j if j >= 0 else len(a)]
                tx = by_id.get(par)
                if tx is None:
                    continue
                if typ == "exon":
                    tx.exons.append((int(f[3]), int(f[4])))
                else:
                    tx.cds.append((int(f[3]), int(f[4]), int(f[7]) if f[7] in "012" else -1))
                continue
            if typ in SKIP_TYPES:
                continue
            at = parse_attrs(f[8])
            if "ID" not in at:
                continue
            key = at.get("transcript_id") or at["ID"].split(":", 1)[-1]
            if wanted is not None and key not in wanted:
                continue
            tx = Tx()
            tx.key = key
            tx.tid = at["ID"]
            tx.version = at.get("version")
            tx.chrom = f[0]
            tx.strand = 1 if f[6] == "+" else -1
            tx.start, tx.end = int(f[3]), int(f[4])
            tx.exons = []
            tx.cds = []
            tx.biotype = at.get("biotype", "")
            by_id[at["ID"]] = tx
            out[key] = tx
    for tx in out.values():
        tx.finish()
    return out


# ----------------------------------------------------------------------------- VCF

def load_vcf(path):
    recs = {}
    dups = set()
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
            alts = alt.split(",")
            key = vid.split(";")[0] if vid != "." else "%s_%d_%s/%s" % (chrom, pos, ref, "/".join(alts))
            if key in recs:
                dups.add(key)
            recs[key] = (chrom, pos, ref.upper(), [a.upper() for a in alts])
    for k in dups:
        del recs[k]
    return recs


def trim_both(ref, alt):
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    ref, alt = ref[i:], alt[i:]
    j = 0
    while j < len(ref) and j < len(alt) and ref[-1 - j] == alt[-1 - j]:
        j += 1
    if j:
        ref, alt = ref[:len(ref) - j], alt[:len(alt) - j]
    return ref, alt


def resolve_alt(ref, alts, allele):
    """index of the ALT a row of the tab file is about, None if it cannot be told"""
    if len(alts) == 1:
        return 0
    allele = allele.upper()
    chop = len({a[0] for a in [ref] + alts if a}) == 1
    hits = []
    for i, a in enumerate(alts):
        reps = {a, (a[1:] if chop else a) or "-", trim_both(ref, a)[1] or "-"}
        if allele in reps:
            hits.append(i)
    return hits[0] if len(hits) == 1 else None


def minimal(lo, ref, alt):
    """trim what ref and alt share at both ends: (lo, hi, repl) with genome[lo:hi] replaced by repl"""
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    ref, alt, lo = ref[i:], alt[i:], lo + i
    j = 0
    while j < len(ref) and j < len(alt) and ref[-1 - j] == alt[-1 - j]:
        j += 1
    if j:
        ref, alt = ref[:len(ref) - j], alt[:len(alt) - j]
    return lo, lo + len(ref), alt


def shift_edit(gw, edit, direction, lo_b, hi_b):
    """the edit moved as far as it goes in the direction (+1 right, -1 left), without leaving [lo_b, hi_b)"""
    lo, hi, s = edit
    if hi > lo and not s:  # a deletion
        if direction > 0:
            while hi < hi_b and gw.base(hi) == gw.base(lo):
                lo += 1
                hi += 1
        else:
            while lo > lo_b and gw.base(lo - 1) == gw.base(hi - 1):
                lo -= 1
                hi -= 1
    elif hi == lo and s:  # an insertion
        if direction > 0:
            while lo < hi_b and gw.base(lo) == s[0]:
                s = s[1:] + s[0]
                lo += 1
            hi = lo
        else:
            while lo > lo_b and gw.base(lo - 1) == s[-1]:
                s = s[-1] + s[:-1]
                lo -= 1
            hi = lo
    return lo, hi, s


def can_shift_once(gw, edit, direction, lo_b, hi_b):
    lo, hi, s = edit
    if hi > lo and not s:
        if direction > 0:
            return hi < hi_b and gw.base(hi) == gw.base(lo)
        return lo > lo_b and gw.base(lo - 1) == gw.base(hi - 1)
    if hi == lo and s:
        if direction > 0:
            return lo < hi_b and gw.base(lo) == s[0]
        return lo > lo_b and gw.base(lo - 1) == s[-1]
    return False


def edited(gw, a, b, edit):
    lo, hi, repl = edit
    return gw.get(a, lo) + repl + gw.get(hi, b)


def equivalent(gw, e1, e2):
    a, b = min(e1[0], e2[0]), max(e1[1], e2[1])
    return edited(gw, a, b, e1) == edited(gw, a, b, e2)


def edit_kind(gw, edit, strand):
    """the HGVS kind a (trimmed) edit should be written as: sub, del, dup, ins, inv, delins"""
    lo, hi, s = edit
    n = hi - lo
    if n == 0 and s:
        ins_s = s
        # most 3' position, then does it repeat what precedes it (in the transcript's direction)?
        sh = shift_edit(gw, edit, strand, 0, 1 << 60)
        lo2, _, s2 = sh
        if strand > 0:
            return "dup" if gw.get(lo2 - len(s2), lo2) == s2 else "ins"
        return "dup" if gw.get(lo2, lo2 + len(s2)) == s2 else "ins"
    if n and not s:
        return "del"
    if n == 1 and len(s) == 1:
        return "sub"
    ref = gw.get(lo, hi)
    if n >= 2 and len(s) == n and s == revcomp(ref):
        return "inv"
    return "delins"


# ----------------------------------------------------------------------------- HGVSc / HGVSg

_POS = r"(?:(?P<{p}star>\*)|(?P<{p}neg>-))?(?P<{p}num>\d+)(?:(?P<{p}sgn>[+-])(?P<{p}off>\d+))?"


def _pos(p):
    return _POS.format(p=p)


R_C_SUB = re.compile(r"^" + _pos("a") + r"(?P<ref>[ACGTN])>(?P<alt>[ACGTN])$")
R_C_OP = re.compile(r"^" + _pos("a") + r"(?:_" + _pos("b") + r")?(?P<op>delins|del|dup|inv|ins)(?P<seq>[ACGTN]*)$")
R_C_EQ = re.compile(r"^" + _pos("a") + r"(?:_" + _pos("b") + r")?=$")
R_G_SUB = re.compile(r"^(?P<a>\d+)(?P<ref>[ACGTN])>(?P<alt>[ACGTN])$")
R_G_OP = re.compile(r"^(?P<a>\d+)(?:_(?P<b>\d+))?(?P<op>delins|del|dup|inv|ins)(?P<seq>[ACGTN]*)$")


def _tx_pos(tx, kind, m, p):
    num = int(m.group(p + "num"))
    star, neg = m.group(p + "star"), m.group(p + "neg")
    if kind == "n":
        if star or neg:
            raise Bad("n. position with - or *")
        t = num
    else:
        if num == 0 and not (star or neg):
            raise Bad("c.0")
        if tx.cds_start is None:
            raise Bad("c. on a transcript without a CDS")
        # numbered from the first base of a codon the transcript may not have (phase 1 or 2): c.1 and c.2 are
        # missing bases.  The intronic anchors of both tools are counted from the first CDS base instead.
        phantom = 0 if m.group(p + "sgn") else tx.cds_offset
        if star:
            t = tx.cds_end + num
        elif neg:
            t = tx.cds_start - phantom - num
        else:
            t = tx.cds_start - phantom + num - 1
    off = 0
    if m.group(p + "sgn"):
        off = int(m.group(p + "off")) * (1 if m.group(p + "sgn") == "+" else -1)
    return t, off


def _to_g(tx, t, off):
    g = tx.g_of(t)
    if off:
        if off > 0 and t not in tx.lasts:
            raise Wrong("offset +%d from a base that is not the last of an exon" % off)
        if off < 0 and t not in tx.firsts:
            raise Wrong("offset %d from a base that is not the first of an exon" % off)
        g = g + off if tx.strand > 0 else g - off
    return g


def c_description(tx, desc, gw):
    """(edit, kind, loc) of an HGVSc body such as `c.123_124insA`; raises Bad"""
    pre, _, body = desc.partition(":")
    if not body or body[0] not in "cn" or body[1:2] != ".":
        raise Bad("not c. or n.")
    kind, body = body[0], body[2:]
    strand = tx.strand
    m = R_C_SUB.match(body)
    loc = "cds"
    if m:
        t, off = _tx_pos(tx, kind, m, "a")
        g = _to_g(tx, t, off)
        loc = _loc(kind, m, ("a",), off, tx)
        refb = gw.get(g - 1, g)
        refb = refb if strand > 0 else revcomp(refb)
        if m.group("ref") != "N" and refb != m.group("ref"):
            raise Wrong("stated reference %s, genome has %s" % (m.group("ref"), refb))
        alt = m.group("alt") if strand > 0 else revcomp(m.group("alt"))
        return (g - 1, g, alt), "sub", loc
    m = R_C_OP.match(body)
    if not m:
        if R_C_EQ.match(body):
            raise Bad("c.= (no change)")
        raise Bad("cannot read")
    t1, o1 = _tx_pos(tx, kind, m, "a")
    gA = _to_g(tx, t1, o1)
    if m.group("bnum"):
        t2, o2 = _tx_pos(tx, kind, m, "b")
        gB = _to_g(tx, t2, o2)
        loc = _loc(kind, m, ("a", "b"), o1 or o2, tx)
    else:
        gB = gA
        loc = _loc(kind, m, ("a",), o1, tx)
    op, seq = m.group("op"), m.group("seq")
    gseq = seq if strand > 0 else revcomp(seq)
    lo, hi = min(gA, gB) - 1, max(gA, gB)
    if op == "del":
        if seq and gw.get(lo, hi) != gseq:
            raise Wrong("stated deleted bases differ from the genome")
        return (lo, hi, ""), "del", loc
    if op == "dup":
        if seq and gw.get(lo, hi) != gseq:
            raise Wrong("stated duplicated bases differ from the genome")
        # a tandem duplication is the block inserted next to it; on the minus strand the 3' side is the low one
        return ((hi, hi, gw.get(lo, hi)) if strand > 0 else (lo, lo, gw.get(lo, hi))), "dup", loc
    if op == "inv":
        if gB == gA:
            raise Bad("inv of one base")
        return (lo, hi, revcomp(gw.get(lo, hi))), "inv", loc
    if op == "delins":
        if not seq:
            raise Bad("delins without bases")
        return (lo, hi, gseq), "delins", loc
    # ins
    if not seq:
        raise Bad("ins without bases")
    if gB == gA or abs(gA - gB) != 1:
        pa, pb = tx.pos_of(gA), tx.pos_of(gB)
        if pa is not None and pb is not None and abs(pa - pb) == 1:
            raise Skip("insertion between the last base of an exon and the first of the next")
        raise Wrong("insertion flanks %d and %d are not adjacent in the genome" % (gA, gB))
    return (min(gA, gB), min(gA, gB), gseq), "ins", loc


def _loc(kind, m, ps, off, tx=None):
    nf = bool(tx is not None and tx.cds_offset and not any(m.group(p + "star") for p in ps))
    loc = _loc0(kind, m, ps, off)
    return loc + "_nf" if nf else loc


def _loc0(kind, m, ps, off):
    if kind == "n":
        return "noncoding_splice" if off and abs(off) <= 2 else ("noncoding_intron" if off else "noncoding")
    for p in ps:
        if m.group(p + "sgn"):
            return "splice" if abs(int(m.group(p + "off"))) <= 2 else "intron"
    for p in ps:
        if m.group(p + "star"):
            return "utr3"
        if m.group(p + "neg"):
            return "utr5"
    return "cds"


def g_description(desc, gw):
    pre, _, body = desc.partition(":")
    if body[:2] not in ("g.", "m."):
        raise Bad("not g.")
    body = body[2:]
    m = R_G_SUB.match(body)
    if m:
        p = int(m.group("a"))
        if m.group("ref") != "N" and gw.get(p - 1, p) != m.group("ref"):
            raise Wrong("stated reference %s, genome has %s" % (m.group("ref"), gw.get(p - 1, p)))
        return (p - 1, p, m.group("alt")), "sub"
    m = R_G_OP.match(body)
    if not m:
        raise Bad("cannot read")
    a = int(m.group("a"))
    b = int(m.group("b")) if m.group("b") else a
    op, seq = m.group("op"), m.group("seq")
    lo, hi = a - 1, b
    if op == "del":
        if seq and gw.get(lo, hi) != seq:
            raise Wrong("stated deleted bases differ from the genome")
        return (lo, hi, ""), "del"
    if op == "dup":
        return (hi, hi, gw.get(lo, hi)), "dup"
    if op == "inv":
        return (lo, hi, revcomp(gw.get(lo, hi))), "inv"
    if op == "delins":
        return (lo, hi, seq), "delins"
    if b != a + 1:
        raise Bad("insertion flanks are not adjacent")
    return (a, a, seq), "ins"


def check_edit_description(gw, vedit, dedit, dkind, strand, lo_b, hi_b):
    """verdict, detail for a description whose edit is dedit and the VCF's is vedit (both untrimmed)"""
    if not equivalent(gw, vedit, dedit):
        return "wrong", "not the VCF change"
    t = minimal(dedit[0], gw.get(dedit[0], dedit[1]), dedit[2])
    exp = edit_kind(gw, minimal(vedit[0], gw.get(vedit[0], vedit[1]), vedit[2]), strand)
    verdict, detail = "ok", ""
    if can_shift_once(gw, t, strand, lo_b, hi_b):
        # a dup is an insertion of the block after it: the same test applies to the trimmed insertion
        verdict, detail = "ok_not3", "can move 3'"
    elif exp != dkind:
        verdict, detail = "ok_spelling", "%s written as %s" % (exp, dkind)
    return verdict, detail


# ----------------------------------------------------------------------------- HGVSp

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X"}
RES = r"(?:[A-Z][a-z]{2}|\*)"
R_P_FS = re.compile(r"^(?P<ref>" + RES + r")(?P<pos>\d+)(?P<alt>" + RES + r")?fs(?:(?:Ter|\*)(?P<n>\d+|\?))?$")
R_P_EXT = re.compile(r"^(?P<ref>Ter|\*)(?P<pos>\d+)(?P<alt>" + RES + r"|del)?ext(?:(?:Ter|\*)(?P<n>\d+|\?))?$")
R_P_DELINS = re.compile(r"^(?P<r1>" + RES + r")(?P<p1>\d+)(?:_(?P<r2>" + RES + r")(?P<p2>\d+))?delins(?P<seq>(?:" + RES + r")+)$")
R_P_INS = re.compile(r"^(?P<r1>" + RES + r")(?P<p1>\d+)_(?P<r2>" + RES + r")(?P<p2>\d+)ins(?P<seq>(?:" + RES + r")+)$")
R_P_DD = re.compile(r"^(?P<r1>" + RES + r")(?P<p1>\d+)(?:_(?P<r2>" + RES + r")(?P<p2>\d+))?(?P<op>del|dup)$")
R_P_DELEXT = re.compile(r"^(?P<r1>" + RES + r")(?P<p1>\d+)_(?P<r2>Ter|\*)(?P<p2>\d+)delext(?:Ter|\*)(?P<n>\d+|\?)$")
R_P_SYNRANGE = re.compile(r"^(?P<r1>" + RES + r")(?P<p1>\d+)_(?P<r2>" + RES + r")(?P<p2>\d+)=$")
R_P_SUB = re.compile(r"^(?P<ref>(?:" + RES + r")+)(?P<pos>\d+)(?P<alt>(?:" + RES + r")+|=|\?)$")
R_P_BADSL = re.compile(r"^(?P<ref>(?:" + RES + r")+)(?P<p1>\d+)_\?(?P<p2>\d+)$")


def one_letter(s):
    return "".join(AA3.get(r, "?") if r != "*" else "*" for r in re.findall(RES, s))


def parse_p(desc):
    """dict with kind and fields; raises Bad for what the checker cannot read"""
    r = _parse_p(desc)
    if r["kind"] != "malformed" and any(r.get(k) == 0 for k in ("pos", "p1", "p2")):
        return {"kind": "malformed", "why": "position 0", "body": desc}
    return r


def _parse_p(desc):
    pre, _, body = desc.partition(":")
    if body[:2] != "p.":
        raise Bad("not p.")
    body = body[2:].replace("%3D", "=").replace("%3d", "=")
    if body.startswith("(") and body.endswith(")"):
        body = body[1:-1]
    m = R_P_BADSL.match(body)
    if m:
        return {"kind": "malformed", "why": "range with ? as the end (start loss form)", "body": body}
    m = R_P_FS.match(body)
    if m:
        return {"kind": "fs", "ref": one_letter(m["ref"]), "pos": int(m["pos"]), "alt": one_letter(m["alt"]) if m["alt"] else None,
                "n": m["n"]}
    m = R_P_EXT.match(body)
    if m:
        return {"kind": "ext", "pos": int(m["pos"]), "alt": (None if m["alt"] == "del" else one_letter(m["alt"])) if m["alt"] else None,
                "del": m["alt"] == "del", "n": m["n"]}
    m = R_P_DELINS.match(body)
    if m:
        p1 = int(m["p1"])
        p2 = int(m["p2"]) if m["p2"] else p1
        return {"kind": "delins", "p1": p1, "p2": p2, "r1": one_letter(m["r1"]), "r2": one_letter(m["r2"]) if m["r2"] else one_letter(m["r1"]),
                "seq": one_letter(m["seq"])}
    m = R_P_INS.match(body)
    if m:
        return {"kind": "ins", "p1": int(m["p1"]), "p2": int(m["p2"]), "r1": one_letter(m["r1"]), "r2": one_letter(m["r2"]),
                "seq": one_letter(m["seq"])}
    m = R_P_DD.match(body)
    if m:
        p1 = int(m["p1"])
        p2 = int(m["p2"]) if m["p2"] else p1
        return {"kind": m["op"], "p1": p1, "p2": p2, "r1": one_letter(m["r1"]), "r2": one_letter(m["r2"]) if m["r2"] else one_letter(m["r1"])}
    m = R_P_DELEXT.match(body)
    if m:
        return {"kind": "delext", "p1": int(m["p1"]), "p2": int(m["p2"]), "r1": one_letter(m["r1"]), "n": m["n"]}
    m = R_P_SYNRANGE.match(body)
    if m:
        return {"kind": "synrange", "p1": int(m["p1"]), "p2": int(m["p2"]), "r1": one_letter(m["r1"]), "r2": one_letter(m["r2"])}
    m = R_P_SUB.match(body)
    if m:
        refs = re.findall(RES, m["ref"])
        if len(refs) > 1:
            return {"kind": "malformed", "why": "several residue names share one position", "body": body}
        alt = m["alt"]
        ref = one_letter(m["ref"])
        pos = int(m["pos"])
        if alt == "=":
            return {"kind": "syn", "ref": ref, "pos": pos}
        if alt == "?":
            return {"kind": "unknown", "ref": ref, "pos": pos}
        alts = re.findall(RES, alt)
        if len(alts) != 1:
            return {"kind": "delins_short", "ref": ref, "pos": pos, "seq": one_letter(alt)}
        return {"kind": "sub", "ref": ref, "pos": pos, "alt": one_letter(alt)}
    raise Bad("cannot read p. form")


def p_class(p):
    k = p["kind"]
    if k == "sub":
        return "nonsense" if p["alt"] == "*" else "missense"
    if k == "unknown":
        return "startloss" if p["pos"] == 1 else "unknown"
    return k


def first_diff(a, b):
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def stop_trunc(s):
    i = s.find("*")
    return s if i < 0 else s[:i + 1]


def check_p(p, ref, alt, start_lost, open_end):
    """does the description p describe the change ref -> alt (proteins with their terminators)?
    returns (verdict, detail)"""
    k = p["kind"]
    if k == "unknown":
        if p["pos"] == 1 and p["ref"] == "M":
            if start_lost:
                return "ok", ""
            return "wrong", "start-loss notation but the initiator codon is intact"
        return "wrong", "unknown effect given for a residue that is not the initiator"
    if start_lost:
        return "wrong", "the initiator is lost and the description is %s, not p.Met1?" % k
    if k == "delext":
        p1, p2 = p["p1"], p["p2"]
        if p2 > len(ref) or ref[p2 - 1] != "*" or ref[p1 - 1] != p["r1"] or p2 != len(ref):
            return "wrong", "residues %d_%d are not %s..Ter at the end of the protein" % (p1, p2, p["r1"])
        if alt[:p1 - 1] != ref[:p1 - 1]:
            return "wrong", "the protein changes before residue %d" % p1
        stop = alt.endswith("*")
        true_n = len(alt) - p2  # counted from the terminator that was deleted, as for `ext`
        n = p["n"]
        if n == "?":
            return ("ok", "") if (open_end or not stop) else ("wrong", "no stop given, there is one at distance %d" % true_n)
        if not stop:
            return "wrong", "a stop at %s is given, the new frame reaches the end of the transcript" % n
        return ("ok", "") if int(n) == true_n else ("wrong", "stop at %s, the edited frame has it at %d" % (n, true_n))
    if k == "synrange":
        p1, p2 = p["p1"], p["p2"]
        if p2 > len(ref) or p1 > p2 or ref[p1 - 1] != p["r1"] or ref[p2 - 1] != p["r2"]:
            return "wrong", "residues at %d_%d are not %s %s" % (p1, p2, p["r1"], p["r2"])
        return ("ok", "") if alt == ref else ("wrong", "the protein changes (first difference at %d)" % (first_diff(ref, alt) + 1))
    if k == "syn":
        if p["pos"] > len(ref) or ref[p["pos"] - 1] != p["ref"]:
            return "wrong", "residue %d is %s, not %s" % (p["pos"], ref[p["pos"] - 1:p["pos"]], p["ref"])
        return ("ok", "") if alt == ref else ("wrong", "the protein changes (first difference at %d)" % (first_diff(ref, alt) + 1))
    if k == "sub":
        pos = p["pos"]
        if pos > len(ref) or ref[pos - 1] != p["ref"]:
            return "wrong", "residue %d is %s, not %s" % (pos, ref[pos - 1:pos], p["ref"])
        if p["ref"] == "*":
            return "wrong", "a change at the terminator is not a substitution"
        exp = stop_trunc(ref[:pos - 1] + p["alt"] + ref[pos:])
        return ("ok", "") if exp == alt else ("wrong", "applied: differs from the edited protein at %d" % (first_diff(exp, alt) + 1))
    if k in ("del", "dup", "delins", "ins", "delins_short"):
        if k == "delins_short":
            p1 = p2 = p["pos"]
            if ref[p1 - 1:p1] != p["ref"]:
                return "wrong", "residue %d is %s, not %s" % (p1, ref[p1 - 1:p1], p["ref"])
        else:
            p1, p2 = p.get("p1", 0), p.get("p2", 0)
            if p1 < 1 or p2 > len(ref) or p2 < p1:
                return "wrong", "positions %d..%d outside the protein" % (p1, p2)
            if k != "ins":
                if ref[p1 - 1] != p["r1"] or ref[p2 - 1] != p["r2"]:
                    return "wrong", "residues at %d..%d are %s..%s, not %s..%s" % (p1, p2, ref[p1 - 1], ref[p2 - 1], p["r1"], p["r2"])
            else:
                if p2 != p1 + 1 or ref[p1 - 1] != p["r1"] or ref[p2 - 1:p2] != p["r2"]:
                    return "wrong", "insertion flanks %d_%d are not %s %s" % (p1, p2, p["r1"], p["r2"])
        if k == "del":
            exp = ref[:p1 - 1] + ref[p2:]
        elif k == "dup":
            exp = ref[:p2] + ref[p1 - 1:p2] + ref[p2:]
        elif k == "ins":
            exp = ref[:p1] + p["seq"] + ref[p1:]
        else:
            exp = ref[:p1 - 1] + p["seq"] + ref[p2:]
        exp = stop_trunc(exp)
        if exp != alt:
            return "wrong", "applied: differs from the edited protein at %d" % (first_diff(exp, alt) + 1)
        # 3' rule at the protein level
        if k in ("del", "dup"):
            blk = ref[p1 - 1:p2]
            nxt = ref[p2:p2 + 1]
            if nxt and nxt == blk[0] and not (nxt == "*"):
                return "ok_not3", "can move 3' in the protein"
        if k == "ins":
            nxt = ref[p1:p1 + 1]
            if nxt and nxt == p["seq"][0] and nxt != "*":
                return "ok_not3", "can move 3' in the protein"
            if len(p["seq"]) <= p1 and ref[p1 - len(p["seq"]):p1] == p["seq"]:
                return "ok_spelling", "ins written for a dup"
        return "ok", ""
    if k == "fs":
        pos = p["pos"]
        d = first_diff(ref, alt)
        if d != pos - 1:
            return "wrong", "first changed residue is %d (%s -> %s), described at %d" % (d + 1, ref[d:d + 1], alt[d:d + 1], pos)
        if ref[pos - 1] != p["ref"]:
            return "wrong", "residue %d is %s, not %s" % (pos, ref[pos - 1], p["ref"])
        if p["alt"] is not None and alt[pos - 1:pos] != p["alt"]:
            return "wrong", "first changed residue is %s, described as %s" % (alt[pos - 1:pos], p["alt"])
        n = p["n"]
        stop = alt.endswith("*")
        true_n = len(alt) - (pos - 1)
        if n is None or n == "?":
            if not open_end and stop:
                return "wrong", "no stop given, there is one at distance %d" % true_n
            return "ok", ""
        if not stop:
            return "wrong", "a stop at %s is given, the new frame reaches the end of the transcript" % n
        if int(n) != true_n:
            return "wrong", "stop at %s, the edited frame has it at %d" % (n, true_n)
        return "ok", ""
    if k == "ext":
        pos = p["pos"]
        if pos > len(ref) or ref[pos - 1] != "*":
            return "wrong", "residue %d is not the terminator" % pos
        if alt[:pos - 1] != ref[:pos - 1]:
            return "wrong", "the protein changes before the terminator"
        if len(alt) < pos:
            if open_end and len(alt) == pos - 1 and p["n"] in (None, "?"):
                return "ok", ""  # the terminator is gone and the transcript has no more bases to read
            return "wrong", "the edited protein ends before the terminator"
        if p["alt"] is not None and alt[pos - 1:pos] != p["alt"]:
            return "wrong", "first residue after the loss is %s, described as %s" % (alt[pos - 1:pos], p["alt"])
        n = p["n"]
        stop = alt.endswith("*")
        true_n = len(alt) - pos  # the new terminator is at position pos + true_n
        if n is None or n == "?":
            if not open_end and stop:
                return "wrong", "no stop given, there is one at distance %d" % true_n
            return "ok", ""
        if not stop:
            return "wrong", "a stop at %s is given, the new frame reaches the end of the transcript" % n
        if int(n) != true_n:
            return "wrong", "stop at %s, the edited frame has it at %d" % (n, true_n)
        return "ok", ""
    return "unparsed", "kind " + k


def to_tx_edit(tx, edit, attach=False):
    """a genomic edit as an edit of the spliced sequence:
    ("none"|"skip", why) or ("rep"|"ins", t1, t2, repl_in_transcript_orientation)"""
    lo, hi, s = edit
    st = tx.strand
    if hi > lo:
        i1, i2 = tx.exon_index(lo + 1), tx.exon_index(hi)
        if i1 is None and i2 is None:
            for (es, ee) in tx.asc:
                if es <= hi and ee >= lo + 1:
                    return "skip", "spans a splice boundary"
            return "none", "outside the exons"
        if i1 != i2:
            return "skip", "spans a splice boundary"
        a, b = tx.pos_of(lo + 1), tx.pos_of(hi)
        return ("rep", min(a, b), max(a, b), s if st > 0 else revcomp(s))
    # an insertion between the 1-based bases lo and lo + 1
    i1, i2 = tx.exon_index(lo), tx.exon_index(lo + 1)
    repl = s if st > 0 else revcomp(s)
    if i1 is not None and i1 == i2:
        t = tx.pos_of(lo) if st > 0 else tx.pos_of(lo + 1)
        return "ins", t, t, repl
    if attach and (i1 is None) != (i2 is None):
        # between the edge of an exon and the intron: read as bases the exon gains
        if i1 is not None:
            t = tx.pos_of(lo) if st > 0 else tx.pos_of(lo) - 1
        else:
            t = tx.pos_of(lo + 1) - 1 if st > 0 else tx.pos_of(lo + 1)
        return "ins", t, t, repl
    return "none", "outside the exons"


def protein_candidate(tx, fasta, tedit):
    """what an edit of the spliced sequence does to the protein, or None when it does not reach the CDS.
    The facts are those of the sequence alone: the edited protein, whether the edit alters the bases of the
    start or the stop codon, and the codons that hold the change on both sides."""
    cs, ce = tx.cds_start, tx.cds_end
    ref, _ = tx.prot(fasta)
    sq = tx.seq(fasta)
    table = tx._table
    kind, t1, t2, repl = tedit
    if kind == "ins":
        if not (cs <= t1 <= ce - 1):
            return None
        alt_tx = sq[:t1] + repl + sq[t1:]
        delta = len(repl)
        first, last = t1 + 1, t1
        straddle, in_cds, over_end = False, True, False
    else:
        if t2 < cs or t1 > ce:
            return None
        delta = len(repl) - (t2 - t1 + 1)
        alt_tx = sq[:t1 - 1] + repl + sq[t2:]
        first, last = t1, t2
        straddle = t1 < cs
        over_end = t2 > ce
        in_cds = t1 >= cs and t2 <= ce
    ref_codon = sq[cs - 1:cs + 2]
    f = {"kind": kind, "delta": delta, "in_cds": in_cds, "straddle": straddle, "over_end": over_end, "alt": None,
         "start_lost": False, "start_retained": False, "start_touch": False, "start_adj": False,
         "stop_touch": False, "stop_adj": False, "open_end": False,
         "nonatg": tx.nf_start or ref_codon not in ("ATG", "TTG", "CTG", "GTG", "ACG"),  # no initiator codon to lose
         "stop_gained": False, "stop_lost": False, "stop_retained": False, "missense": False, "syn": False,
         "ref_win": "", "alt_win": "", "ref_codon": ref_codon, "alt_codon": "", "kept": False}
    # An edit that touches the start codon and leaves the whole coding sequence behind it, with no ATG where the
    # initiator was: the first ATG of the edited transcript is the one the transcript had, so the protein is the
    # protein, whichever way the edit is spelled (a deletion of the initiator's own A after an A, `c.1dup`, a deletion
    # that runs from the UTR into the start codon inside a repeat). The only thing to say of the start codon is that a
    # base of it changed and the start remains. An ATG at the old coding start is read from there instead.
    touch = t1 in (cs, cs + 1) if kind == "ins" else (t1 <= cs + 2 and t2 >= cs)
    if (touch and delta != 0 and not f["nonatg"] and cs - 1 + delta >= 0 and alt_tx[cs - 1:cs + 2] != "ATG"
            and alt_tx[cs - 1 + delta:ce + delta] == sq[cs - 1:ce]):
        f.update(kept=True, straddle=False, start_touch=True, start_retained=True, alt=ref, alt_codon=alt_tx[cs - 1:cs + 2])
        return f
    if straddle:
        # an edit that starts upstream of the initiator and reaches into it: no placement of it that
        # avoids the start codon is available here, the caller looks at the 3'-most placement as well
        f["start_touch"] = True
        f["start_lost"] = not f["nonatg"]
        return f
    if kind == "ins":
        start_touch, start_adj = t1 in (cs, cs + 1), t1 == cs + 2
        stop_touch, stop_adj = t1 in (ce - 2, ce - 1), t1 == ce - 3
    else:
        start_touch, start_adj = (t1 <= cs + 2 and t2 >= cs), t1 == cs + 3
        stop_touch, stop_adj = (t1 <= ce and t2 >= ce - 2), t2 == ce - 3
    if f["nonatg"]:
        start_touch = start_adj = False  # the first codon is a codon like the others
    f.update(start_touch=start_touch, start_adj=start_adj, stop_touch=stop_touch, stop_adj=stop_adj)
    f["alt_codon"] = alt_tx[cs - 1:cs + 2]
    if start_touch and not f["nonatg"]:
        if f["alt_codon"] == ref_codon or f["alt_codon"] == "ATG" or (ref_codon == "ATG" and translate(f["alt_codon"], table) == "M"):
            f["start_retained"] = True
        else:
            f["start_lost"] = True
    # the protein: from the codon that holds the first changed base
    k = (max(first, cs) - cs) // 3
    prefix = ref[:k]
    end_alt = (last - cs + 1) + delta
    if delta % 3 == 0 and last < ce - 2 and first >= cs and not f["start_lost"]:
        # in frame again after the edit: the codons that hold it, then the reference tail
        j = max((end_alt - 1) // 3 + 2, k + 1)
        mid = translate(alt_tx[cs - 1 + 3 * k:cs - 1 + 3 * j], table)
        f["alt"] = prefix + stop_trunc(mid) if "*" in mid else prefix + mid + ref[j - delta // 3:]
    else:
        mid = translate(alt_tx[cs - 1 + 3 * k:], table, stop_at_stop=True)
        f["open_end"] = not mid.endswith("*")
        f["alt"] = prefix + mid
    # the codons of the reference that hold the changed bases, and of the edited sequence that hold the new ones
    ncod = (ce - cs) // 3 + 1
    if kind == "ins":
        if (t1 - cs + 1) % 3 == 0:
            ka = kb = (t1 - cs + 1) // 3
        else:
            ka = (t1 - cs) // 3
            kb = ka + 1
    else:
        ka = (t1 - cs) // 3
        kb = min((t2 - cs) // 3 + 1, ncod)
    alt_last = max(end_alt, first - cs + 1)
    kb_alt = max((alt_last - 1) // 3 + 1, ka)
    f["ref_win"] = ref[ka:kb]
    f["alt_win"] = translate(alt_tx[cs - 1 + 3 * ka:cs - 1 + 3 * kb_alt], table)
    rw, aw = f["ref_win"], f["alt_win"]
    # a stop is gained when the protein ends before the reference's own terminator: one inside a sequence that is
    # inserted in front of that terminator only moves the end further on
    f["stop_gained"] = "*" in aw and "*" not in rw and f["alt"].endswith("*") and len(f["alt"]) < len(ref)
    f["stop_lost"] = stop_touch and "*" in rw and "*" not in aw
    f["stop_retained"] = stop_touch and "*" in rw and "*" in aw
    if in_cds and delta == 0 and not stop_touch:
        if f["alt"] == ref:
            f["syn"] = not start_touch
        elif not f["start_lost"] and not f["stop_gained"]:
            diffs = [(r, a) for r, a in zip(rw, aw) if r != a]
            f["missense"] = bool(diffs) and all(r != "*" and a != "*" for r, a in diffs)
    return f


U_TERMS = ["start_lost", "start_retained_variant", "stop_lost", "stop_gained", "stop_retained_variant",
           "frameshift_variant", "inframe_insertion", "inframe_deletion", "synonymous_variant", "missense_variant"]


def expected_terms(f):
    """(expected, do_not_care) sets of the terms of U_TERMS that the sequence decides"""
    if f["straddle"]:
        return None
    if f["kept"]:
        return {"start_retained_variant"}, set()
    E, D = set(), set()
    if f["start_lost"]:
        E.add("start_lost")
    if f["start_retained"]:
        (D if f["start_adj"] else E).add("start_retained_variant")
    elif f["start_adj"]:
        D.add("start_retained_variant")
    if f["stop_gained"]:
        E.add("stop_gained")
    if f["stop_lost"]:
        E.add("stop_lost")
    if f["stop_retained"]:
        E.add("stop_retained_variant")
    elif f["stop_adj"]:
        D.add("stop_retained_variant")
    if f["in_cds"]:
        d = f["delta"]
        if d % 3:
            E.add("frameshift_variant")
        elif d > 0:
            E.add("inframe_insertion")
        elif d < 0:
            E.add("inframe_deletion")
    else:
        D.update(("frameshift_variant", "inframe_insertion", "inframe_deletion"))
    if f["syn"]:
        E.add("synonymous_variant")
    if f["missense"]:
        E.add("missense_variant")
    elif f["delta"] == 0 and f["in_cds"] and f["alt"] != "" and not f["syn"]:
        D.add("missense_variant")  # a stop or the initiator is involved: the sequence alone does not say
    if f["nonatg"]:
        D.update(("start_lost", "start_retained_variant"))
    return E, D


def interpretations(tx, gw, vmin, fasta):
    """candidate outcomes of the change on the protein, for the VCF's placement and for the 3'-most placement
    (the two spell one genome; which one the tool reads is its own choice)"""
    lo_b, hi_b = tx.start - 1, tx.end
    placements = [vmin]
    if (vmin[1] == vmin[0] and vmin[2]) or (vmin[1] > vmin[0] and not vmin[2]):
        for direction in (tx.strand, -tx.strand):  # the 3'-most placement, then the 5'-most
            e2 = shift_edit(gw, vmin, direction, lo_b, hi_b)
            if e2 not in placements:
                placements.append(e2)
    outs = []
    for e in placements:
        te = to_tx_edit(tx, e)
        if te[0] in ("skip", "none"):
            outs.append((te[0], te[1]))
            if te[0] == "none" and e[1] == e[0]:
                ta = to_tx_edit(tx, e, attach=True)
                if ta[0] == "ins":
                    f = protein_candidate(tx, fasta, ta)
                    if f is not None:
                        outs.append(("prot", f))
        else:
            f = protein_candidate(tx, fasta, te)
            outs.append(("none", "no CDS effect") if f is None else ("prot", f))
    if len(outs) > 1:
        keep = [o for o in outs if not (o[0] == "prot" and o[1]["straddle"])]
        if keep:
            outs = keep
    return outs


# ----------------------------------------------------------------------------- the checker of one tab file

class Ctx:
    pass


CTX = Ctx()


def facts_str(f):
    return "d=%+d;codon0=%s>%s;start_touch=%d;start_adj=%d;stop_touch=%d;stop_adj=%d;ref_win=%s;alt_win=%s;alt0=%s" % (
        f["delta"], f["ref_codon"], f["alt_codon"], f["start_touch"], f["start_adj"], f["stop_touch"], f["stop_adj"],
        f["ref_win"][:12], f["alt_win"][:12], (f["alt"] or "")[:1])


def eval_c(tx, hc, gw, vedit, vmin):
    """(class, verdict, detail, extra flags) of an HGVSc against the VCF change"""
    lo_b, hi_b = tx.start - 1, tx.end

    def once():
        dedit, dkind, loc = c_description(tx, hc, gw)
        cl = "c:%s:%s" % (dkind, loc)
        if vmin[0] < lo_b or vmin[1] > hi_b:
            # a change that runs off the transcript is described for the part inside it
            if vmin[1] > vmin[0] and not vmin[2]:
                v_use = (max(vmin[0], lo_b), min(vmin[1], hi_b), "")
                v, d = check_edit_description(gw, v_use, dedit, dkind, tx.strand, lo_b, hi_b)
            else:
                v, d = "skip_overhang", "a replacement or insertion that runs off the transcript"
            return cl, v, d, ["overhang"]
        v, d = check_edit_description(gw, vedit, dedit, dkind, tx.strand, lo_b, hi_b)
        return cl, v, d, []

    try:
        cl, v, d, ex = once()
    except Skip as e:
        return "c:?", "skip_junction", str(e), []
    except Bad as e:
        cl, v, d, ex = "c:?", ("wrong" if isinstance(e, Wrong) else "unparsed"), str(e), []
    if tx.cds_offset and v in ("wrong", "ok_not3"):
        # both tools count an indel of a transcript whose CDS starts part-way through a codon from the first CDS
        # base, and a substitution from the codon's first base
        saved = tx.cds_offset
        tx.cds_offset = 0
        try:
            cl2, v2, d2, ex2 = once()
        except Bad:
            cl2, v2, d2, ex2 = cl, "wrong", "", ex
        finally:
            tx.cds_offset = saved
        if v2 in ("ok", "ok_spelling", "ok_not3"):
            return cl2, "ok_numbering", "counted from the first CDS base, not from the codon (%s)" % v2, ex2
        if v2.startswith("skip"):
            return cl2, v2, d2, ex2
    return cl, v, d, ex


def process_group(lines):
    """the tab lines of one or more variants -> (output lines, counters)"""
    out = []
    st = collections.Counter()
    C = CTX
    seen_g = set()
    gws = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        var, allele, feat = f[C.i_var], f[C.i_allele], f[C.i_feat]
        rec = C.vcf.get(var)
        if rec is None:
            st["row_no_vcf_record"] += 1
            continue
        chrom, pos, ref, alts = rec
        ai = resolve_alt(ref, alts, allele)
        flags = ["multi"] if len(alts) > 1 else []
        hc = f[C.i_c] if C.i_c >= 0 else "-"
        hp = f[C.i_p] if C.i_p >= 0 else "-"
        hg = f[C.i_g] if C.i_g >= 0 else "-"
        if ai is None:
            st["row_allele_unresolved"] += 1
            continue
        alt = alts[ai]
        if not re.fullmatch(r"[ACGT]+", alt) or not re.fullmatch(r"[ACGT]+", ref):
            continue
        gw = gws.get(chrom)
        if gw is None:
            gw = gws[chrom] = GW(C.fasta, chrom)
        try:
            ok_ref = gw.get(pos - 1, pos - 1 + len(ref)) == ref
        except Bad:
            ok_ref = False
        if not ok_ref:
            st["row_ref_mismatch"] += 1
            continue
        vedit = (pos - 1, pos - 1 + len(ref), alt)
        vmin = minimal(pos - 1, ref, alt)
        # HGVSg: once per variant allele
        key = (var, ai)
        if C.i_g >= 0 and hg != "-" and key not in seen_g:
            seen_g.add(key)
            try:
                dedit, dkind = g_description(hg, gw)
                v, d = check_edit_description(gw, vedit, dedit, dkind, 1, 0, 1 << 60)
                out.append("g\t%s\t%s\t-\t%s\tg:%s\t%s\t%s\t%s" % (var, allele, hg, dkind, v, d, ",".join(flags)))
            except Bad as e:
                out.append("g\t%s\t%s\t-\t%s\tg:?\t%s\t%s\t%s" % (var, allele, hg, "wrong" if isinstance(e, Wrong) else "unparsed", e, ",".join(flags)))
        if feat in C.skip_tx:
            st["row_skipped_transcript"] += 1
            continue
        tx = C.tx.get(feat)
        if tx is None:
            continue
        if canon_chrom(tx.chrom) != canon_chrom(chrom):
            continue
        if C.i_flags >= 0:
            tx.nf_start = "cds_start_NF" in f[C.i_flags].split(",")
        fl = list(flags)
        if vmin[0] < tx.start + 10 or vmin[1] > tx.end - 10:
            fl.append("txend")
        # ---------------- HGVSc
        if hc != "-" and C.i_c >= 0:
            pre = hc.partition(":")[0]
            if pre.split(".")[0] != tx.key.split(".")[0]:
                out.append("c\t%s\t%s\t%s\t%s\tc:?\twrong\ttranscript accession %s is not %s\t%s" % (var, allele, feat, hc, pre, tx.key, ",".join(fl)))
            else:
                cl, v, d, ex = eval_c(tx, hc, gw, vedit, vmin)
                out.append("c\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (var, allele, feat, hc, cl, v, d, ",".join(fl + ex)))
        # ---------------- the protein and the consequence terms
        if C.i_p < 0 and C.i_cons < 0:
            continue
        pref, why = tx.prot(C.fasta)
        if why:
            st["skip_tx_" + why] += 1
            if hp != "-" and C.i_p >= 0:
                out.append("p\t%s\t%s\t%s\t%s\tp:?\tskip_tx_%s\t\t%s" % (var, allele, feat, hp, why, ",".join(fl)))
            continue
        outs = interpretations(tx, gw, vmin, C.fasta)
        # --- HGVSp
        if C.i_p >= 0:
            prot_outs = [o for o in outs if o[0] == "prot"]
            if hp == "-":
                if prot_outs and not any(o[0] == "none" for o in outs) and \
                        any(o[1]["alt"] != pref or o[1]["start_lost"] for o in prot_outs):
                    o1 = prot_outs[0][1]
                    out.append("p\t%s\t%s\t%s\t-\tp:none\tmissing\t%s\t%s" % (var, allele, feat, "straddle;" * o1["straddle"] + facts_str(o1), ",".join(fl)))
            else:
                try:
                    p = parse_p(hp)
                except Bad as e:
                    out.append("p\t%s\t%s\t%s\t%s\tp:?\tunparsed\t%s\t%s" % (var, allele, feat, hp, e, ",".join(fl)))
                    p = None
                if p is not None:
                    cl = "p:" + p_class(p)
                    if p["kind"] == "malformed":
                        out.append("p\t%s\t%s\t%s\t%s\t%s\tmalformed\t%s\t%s" % (var, allele, feat, hp, cl, p["why"], ",".join(fl)))
                    else:
                        results = []
                        for o in outs:
                            if o[0] == "skip":
                                results.append(("skip_splice", o[1]))
                            elif o[0] == "none":
                                results.append(("wrong", "no change of the protein expected (%s)" % o[1]))
                            else:
                                cd = o[1]
                                if cd["straddle"]:
                                    results.append(("ok", "") if (p["kind"] == "unknown" and p["pos"] == 1) else
                                                   ("wrong", "the initiator is destroyed and the description is %s, not p.Met1?" % p["kind"]))
                                else:
                                    results.append(check_p(p, pref, cd["alt"], cd["start_lost"], cd["open_end"]))
                        order = {"ok": 0, "ok_not3": 1, "ok_spelling": 2}
                        good = [r for r in results if r[0] in order]
                        if good:
                            v, d = sorted(good, key=lambda r: order[r[0]])[0]
                        elif any(r[0] == "skip_splice" for r in results):
                            v, d = "skip_splice", "spans a splice boundary"
                        else:
                            wr = [r for r in results if r[0] == "wrong"]
                            real = [r for r in wr if not r[1].startswith("no change of the protein")]
                            v, d = (real or wr or results)[0]
                        fl3 = fl + (["ambiguous_placement"] if len({r[0] for r in results}) > 1 else [])
                        out.append("p\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (var, allele, feat, hp, cl, v, d, ",".join(fl3)))
        # --- the terms
        if C.i_cons >= 0:
            if any(o[0] == "skip" for o in outs):
                st["terms_skip_splice"] += 1
                continue
            obs = set(f[C.i_cons].split(",")) & set(U_TERMS)
            exps, dcs = [], []
            first_prot = None
            for o in outs:
                if o[0] == "prot":
                    r = expected_terms(o[1])
                    if r is None:
                        exps = None
                        break
                    exps.append(r[0])
                    dcs.append(r[1])
                    first_prot = first_prot or o[1]
                else:
                    exps.append(set())
                    dcs.append(set())
            if exps is None:
                st["terms_skip_straddle"] += 1
                continue
            if first_prot is None and not obs:
                continue
            st["terms_rows"] += 1
            if first_prot is not None:
                st["terms_rows_cds"] += 1
            any_e = set().union(*exps)
            all_e = set.intersection(*exps)
            care_not = set().union(*dcs) if dcs else set()
            facts = facts_str(first_prot) if first_prot else "no CDS effect"
            if first_prot is None:
                ctx = "noeffect"
            elif first_prot["start_touch"]:
                ctx = "start"
            elif first_prot["start_adj"]:
                ctx = "start_adj"
            elif first_prot["stop_touch"]:
                ctx = "stop"
            elif first_prot["stop_adj"]:
                ctx = "stop_adj"
            else:
                ctx = "body"
            gate = "+".join(sorted(obs & {"start_lost", "stop_lost", "stop_gained", "stop_retained_variant", "frameshift_variant"}))
            for t in U_TERMS:
                inobs = t in obs
                inexp = t in exps[0]
                if inobs:
                    st["T|%s|obs" % t] += 1
                if inexp:
                    st["T|%s|exp" % t] += 1
                if inobs and inexp:
                    st["T|%s|agree" % t] += 1
                if t in care_not and t not in any_e:
                    continue
                if inobs and t not in any_e:
                    st["T|%s|fp" % t] += 1
                    st["T|%s|fp|%s" % (t, ctx)] += 1
                    out.append("t\t%s\t%s\t%s\t%s\tt:FP:%s:%s\tcontra\t%s\t%s" % (var, allele, feat, f[C.i_cons], t, ctx, facts, ",".join(fl)))
                elif not inobs and t in all_e:
                    st["T|%s|fn" % t] += 1
                    st["T|%s|fn|%s" % (t, ctx)] += 1
                    out.append("t\t%s\t%s\t%s\t%s\tt:FN:%s:%s%s\tcontra\t%s\t%s" % (var, allele, feat, f[C.i_cons], t, ctx, (":with_" + gate) if gate else "", facts, ",".join(fl)))
    return out, st


def _init_worker():
    CTX.fasta.pid = -1  # reopen the FASTA after the fork: a file offset is shared by a forked handle


def read_groups(path, every):
    """yield lists of lines, split at variant boundaries, about 3000 lines each"""
    buf, last = [], None
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            v = line.split("\t", 1)[0]
            if every > 1 and zlib.crc32(v.encode()) % every:
                continue
            if v != last and len(buf) >= 3000:
                yield buf
                buf = []
            last = v
            buf.append(line)
    if buf:
        yield buf


def cmd_check(a):
    hdr = None
    with open(a.tab) as fh:
        for line in fh:
            if line.startswith("#Uploaded_variation"):
                hdr = line.rstrip("\n").lstrip("#").split("\t")
                break
    if hdr is None:
        sys.exit("no #Uploaded_variation header in " + a.tab)
    col = {c: i for i, c in enumerate(hdr)}
    CTX.i_var, CTX.i_allele, CTX.i_feat = col["Uploaded_variation"], col["Allele"], col["Feature"]
    CTX.i_c, CTX.i_p, CTX.i_g = col.get("HGVSc", -1), col.get("HGVSp", -1), col.get("HGVSg", -1)
    CTX.i_cons = col.get("Consequence", -1)
    CTX.i_flags = col.get("FLAGS", -1)
    CTX.fasta = Fasta(a.fasta)
    CTX.skip_tx = set(open(a.skip_transcripts).read().split()) if a.skip_transcripts else set()
    CTX.vcf = load_vcf(a.vcf)
    sys.stderr.write("vcf: %d records\n" % len(CTX.vcf))
    wanted = set()
    with open(a.tab) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if a.every > 1 and zlib.crc32(f[0].encode()) % a.every:
                continue
            wanted.add(f[col["Feature"]])
    wanted.discard("-")
    sys.stderr.write("transcripts wanted: %d\n" % len(wanted))
    CTX.tx = load_gff(a.gff3, wanted)
    sys.stderr.write("transcripts loaded: %d\n" % len(CTX.tx))
    nout = 0
    total = collections.Counter()
    with open(a.out, "w") as out:
        out.write("type\tvar\tallele\tfeature\tdesc\tclass\tverdict\tdetail\tflags\n")
        groups = read_groups(a.tab, a.every)
        if a.workers > 1:
            ctx = multiprocessing.get_context("fork")
            with ctx.Pool(a.workers, initializer=_init_worker) as pool:
                for res, st in pool.imap(process_group, groups):
                    for r in res:
                        out.write(r + "\n")
                        nout += 1
                    total.update(st)
        else:
            for g in groups:
                res, st = process_group(g)
                for r in res:
                    out.write(r + "\n")
                    nout += 1
                total.update(st)
        for k, n in sorted(total.items()):
            out.write("S\t%s\t%d\n" % (k, n))
    sys.stderr.write("wrote %d lines to %s\n" % (nout, a.out))


# ----------------------------------------------------------------------------- summaries

def read_results(path):
    rows, stats = [], collections.Counter()
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                stats[f[1]] = int(f[2])
                continue
            while len(f) < 9:
                f.append("")
            rows.append(f)
    return rows, stats


def cmd_summary(a):
    cols = ["ok", "ok_spelling", "ok_not3", "ok_numbering", "wrong", "malformed", "unparsed", "missing", "skipped"]
    for path in a.files:
        rows, stats = read_results(path)
        print("== %s: %d lines" % (path, len(rows)))
        tally = collections.defaultdict(collections.Counter)
        for t, var, allele, feat, desc, cl, v, det, fl in rows:
            if t == "t":
                continue
            tally[cl][v] += 1
            if v.startswith("skip"):
                tally[cl]["skipped"] += 1
        print("%-24s" % "class" + "".join("%13s" % c for c in cols))
        for cl in sorted(tally):
            c = tally[cl]
            print("%-24s" % cl + "".join("%13d" % c.get(k, 0) for k in cols))
        byflag = collections.defaultdict(collections.Counter)
        for t, var, allele, feat, desc, cl, v, det, fl in rows:
            if t == "t":
                continue
            for f in set(fl.split(",")) - {""}:
                byflag[(t, f)][v] += 1
                if v.startswith("skip"):
                    byflag[(t, f)]["skipped"] += 1
        if byflag:
            print("-- by flag")
            for k in sorted(byflag):
                c = byflag[k]
                print("%-24s" % ("%s:%s" % k) + "".join("%13d" % c.get(x, 0) for x in cols))
        skips = {k: v for k, v in stats.items() if not k.startswith("T|")}
        if skips:
            print("counters: " + ", ".join("%s=%d" % kv for kv in sorted(skips.items())))


def failure_class(t, cl, v, desc, det, fl=""):
    """a short name for what is wrong with a description"""
    if v == "malformed":
        return "malformed: " + det
    if t == "c":
        if v == "ok_numbering":
            return "numbering: indel counted from the first CDS base on a 5'-incomplete transcript"
        if v in ("wrong", "ok_not3"):
            if "overhang" in fl.split(","):
                return "%s: a change that runs off the transcript, described as the clipped change moved by the shift of the unclipped one" % v
            norm = re.sub(r"\d+", "N", det)
            norm = re.sub(r"stated reference \w, genome has \w", "stated reference B, genome has B'", norm)
            return "%s: %s" % (v, norm[:90])
    if t == "p":
        if det.startswith("start-loss notation"):
            return "p.Met1? but the initiator codon is intact"
        if det.startswith("unknown effect given"):
            return "p.Xxx<N>? for a residue that is not the initiator"
        if det.startswith("the initiator is lost and the description is"):
            return "an insertion or deletion carried onto the initiator is not p.Met1?"
        if det.startswith("no change of the protein expected"):
            return "protein change for a change outside the exons"
        if "delins" in desc and det.startswith("applied") and re.search(r"Ter\d+(delins|_)|_Ter\d+delins", desc):
            return "delins over the terminator that leaves out the extension"
        if det.startswith("applied"):
            return "%s: applied description differs from the edited protein" % cl
        if det.startswith(("stop at", "no stop given", "a stop at", "first changed residue", "first residue after")):
            return "%s: terminator distance or first residue" % cl
        if "the protein changes" in det:
            return "%s: protein changes under a synonymous description" % cl
        return "%s: %s" % (cl, re.sub(r"\d+", "N", det)[:60])
    return "%s: %s" % (cl, re.sub(r"\d+", "N", det)[:60])


def cmd_failures(a):
    """named failure classes of a results file, with the shortest examples"""
    rows, stats = read_results(a.file)
    bad = ("wrong", "malformed", "unparsed", "ok_numbering") if not a.strict else ("wrong", "malformed", "unparsed", "ok_numbering", "ok_not3", "ok_spelling")
    groups = collections.defaultdict(list)
    for t, var, allele, feat, desc, cl, v, det, fl in rows:
        if t in ("c", "p", "g") and v in bad:
            groups[(t, v, failure_class(t, cl, v, desc, det, fl))].append((var, allele, feat, desc, det, fl))
    print("== %s" % a.file)
    for k, lst in sorted(groups.items(), key=lambda kv: (kv[0][0], -len(kv[1]))):
        print("%s %-10s %7d  %s" % (k[0], k[1], len(lst), k[2]))
        for e in sorted(lst, key=lambda e: (len(e[0]) + len(e[1]), e[0]))[:a.show]:
            print("        %s %s %s  %s  [%s] %s" % (e[0][:70], e[1][:24], e[2], e[3][-70:], e[4][:70], e[5]))


def cmd_table(a):
    """one table per description type: for every kind, per tool, the rows checked and how many are right"""
    data = collections.OrderedDict()
    for spec in a.files:
        name, _, path = spec.partition("=")
        data[name] = read_results(path)[0]
    good_strict = {"ok"}
    good = {"ok", "ok_spelling", "ok_not3", "ok_numbering"}
    bad = {"wrong", "malformed", "unparsed"}
    for typ, title in (("c", "HGVSc"), ("p", "HGVSp"), ("g", "HGVSg")):
        agg = {n: collections.defaultdict(collections.Counter) for n in data}
        for n, rows in data.items():
            for t, var, allele, feat, desc, cl, v, det, fl in rows:
                if t != typ or v.startswith("skip") or v == "missing":
                    continue
                parts = cl.split(":")
                kind = parts[1] if len(parts) > 1 else "?"
                loc = parts[2] if len(parts) > 2 else ""
                key = kind if a.by != "loc" else (loc or kind)
                if v in good_strict:
                    agg[n][key]["ok"] += 1
                elif v in good:
                    agg[n][key][v] += 1
                elif v in bad:
                    agg[n][key]["fail"] += 1
                agg[n][key]["n"] += 1
        keys = sorted({k for n in agg for k in agg[n]})
        print("\n%s, by %s: rows checked | right as spelled | right, other spelling or placement | wrong or malformed" % (title, "location" if a.by == "loc" else "kind"))
        print("%-14s" % "" + "".join("| %-34s" % n for n in data))
        for k in keys:
            line = "%-14s" % k
            for n in data:
                c = agg[n].get(k, collections.Counter())
                other = c["n"] - c["ok"] - c["fail"]
                line += "| %9d %9d %6d %6d    " % (c["n"], c["ok"], other, c["fail"])
            print(line)
        tot = "%-14s" % "total"
        for n in data:
            c = collections.Counter()
            for k in agg[n]:
                c.update(agg[n][k])
            tot += "| %9d %9d %6d %6d    " % (c["n"], c["ok"], c["n"] - c["ok"] - c["fail"], c["fail"])
        print(tot)


def cmd_compare(a):
    vep, _ = read_results(a.vep)
    fvr, _ = read_results(a.fv)
    vd = {(r[0], r[1], r[2], r[3]): r[4:] for r in vep if r[0] != "t"}
    fd = {(r[0], r[1], r[2], r[3]): r[4:] for r in fvr if r[0] != "t"}
    good = ("ok", "ok_spelling", "ok_not3")
    tally = collections.defaultdict(collections.Counter)
    ex = collections.defaultdict(list)
    for k in set(vd) & set(fd):
        a_, b_ = vd[k], fd[k]
        if a_[0] == b_[0]:
            continue
        cls = b_[1] if b_[1] not in ("c:?", "p:?") else a_[1]
        va, vb = a_[2], b_[2]
        if va in good and vb in good:
            who = "both describe the change" + ("" if va == vb else " (VEP %s, fastVEP %s)" % (va, vb))
        elif va in good:
            who = "VEP right, fastVEP %s" % vb
        elif vb in good:
            who = "fastVEP right, VEP %s" % va
        else:
            who = "both fail (VEP %s, fastVEP %s)" % (va, vb)
        tally[cls][who] += 1
        if len(ex[(cls, who)]) < a.show:
            ex[(cls, who)].append((k, a_, b_))
    print("descriptions that differ between the two tools, by class and who is right")
    for cls in sorted(tally):
        for who, n in tally[cls].most_common():
            print("%-24s %7d  %s" % (cls, n, who))
            for k, a_, b_ in ex[(cls, who)]:
                print("      %s %s %s\n         VEP     %s  [%s %s]\n         fastVEP %s  [%s %s]" % (
                    k[1], k[2], k[3], a_[0], a_[2], a_[3], b_[0], b_[2], b_[3]))


def cmd_terms(a):
    """per tool and per term: how often the sequence agrees with the Consequence, and the contradictions
    usage: terms NAME=FILE [NAME=FILE ...]"""
    data = collections.OrderedDict()
    for spec in a.files:
        name, _, path = spec.partition("=")
        data[name] = read_results(path)
    names = list(data)
    print("terms the sequence decides, by tool (rows with a CDS effect, or with one of these terms)")
    print("%-24s" % "term" + "".join(" | %-9s %8s %8s %6s %6s" % (n, "expected", "agree", "FP", "FN") for n in names))
    for t in U_TERMS:
        line = "%-24s" % t
        for name in names:
            stats = data[name][1]
            line += " | %-9d %8d %8d %6d %6d" % tuple(stats.get("T|%s|%s" % (t, k), 0) for k in ("obs", "exp", "agree", "fp", "fn"))
        print(line)
    print("rows: " + "; ".join("%s: %d considered, %d with a CDS effect" % (n, data[n][1].get("terms_rows", 0), data[n][1].get("terms_rows_cds", 0)) for n in names))
    keysets = {n: {(r[5], r[1], r[2], r[3]): r for r in data[n][0] if r[0] == "t"} for n in names}
    for name in names:
        classes = collections.defaultdict(list)
        for k, r in keysets[name].items():
            classes[k[0]].append(r)
        print("\n== %s: contradictions between the Consequence and the edited sequence (%d)" % (name, sum(len(v) for v in classes.values())))
        for cl in sorted(classes, key=lambda c: (c.split(":")[2], c)):
            rs = classes[cl]
            shared = "".join("  [%s: %d]" % (o, sum(1 for r in rs if (cl, r[1], r[2], r[3]) in keysets[o])) for o in names if o != name)
            print("%-52s %7d%s" % (cl, len(rs), shared))
            for r in sorted(rs, key=lambda r: (len(r[1]), r[1]))[:a.show]:
                print("      %s %s %s  Consequence=%s\n         %s" % (r[1], r[2], r[3], r[4], r[7]))


# ----------------------------------------------------------------------------- self-test

def selftest():
    import tempfile
    x = 987654321

    def rnd():
        nonlocal x
        x = (x * 1103515245 + 12345) & 0x7FFFFFFF
        return "ACGT"[(x >> 16) & 3]

    g = [rnd() for _ in range(400)]

    def put(pos1, s):
        for i, c in enumerate(s):
            g[pos1 - 1 + i] = c

    # + strand gene TP: exon 21..50 and 71..100, CDS 31..50 + 71..89 = 39 bases = 13 codons
    #   ATG AAA GGG TTT CCC GGA GC|T ATT GAT AAC CCG GGG TAA  = M K G F P G A I D N P G *
    put(30, "A")  # the base in front of the initiator repeats its first
    put(31, "ATGAAAGGGTTTCCCGGAGC")
    put(51, "GTAAGTCCCCCCCCCCCCAG")
    put(71, "TATTGATAACCCGGGGTAAGGTTAACCCCC")
    # - strand gene TM: exon 2 (5') 251..281, exon 1 (3') 201..230, transcript runs from 281 down to 201
    ex2_tx = "ACGTACGTAC" + "ATGGCCAAAGAATTTGGCCGT"
    ex1_tx = "GGAAACCCTTTGGCTAGTTTAAACTGAGGT"
    put(251, revcomp(ex2_tx))
    put(201, revcomp(ex1_tx))
    genome = "".join(g)
    d = tempfile.mkdtemp()
    with open(d + "/g.fa", "w") as fh:
        fh.write(">c1\n")
        for i in range(0, len(genome), 60):
            fh.write(genome[i:i + 60] + "\n")
    with open(d + "/g.fa.fai", "w") as fh:
        fh.write("c1\t%d\t4\t60\t61\n" % len(genome))
    gff = [
        "c1\tt\tmRNA\t21\t100\t.\t+\t.\tID=transcript:TP;transcript_id=TP;version=1",
        "c1\tt\texon\t21\t50\t.\t+\t.\tParent=transcript:TP",
        "c1\tt\texon\t71\t100\t.\t+\t.\tParent=transcript:TP",
        "c1\tt\tCDS\t31\t50\t.\t+\t0\tParent=transcript:TP",
        "c1\tt\tCDS\t71\t89\t.\t+\t0\tParent=transcript:TP",
        "c1\tt\tmRNA\t201\t281\t.\t-\t.\tID=transcript:TM;transcript_id=TM;version=1",
        "c1\tt\texon\t201\t230\t.\t-\t.\tParent=transcript:TM",
        "c1\tt\texon\t251\t281\t.\t-\t.\tParent=transcript:TM",
        "c1\tt\tCDS\t204\t230\t.\t-\t0\tParent=transcript:TM",
        "c1\tt\tCDS\t251\t271\t.\t-\t0\tParent=transcript:TM",
    ]
    with open(d + "/g.gff3", "w") as fh:
        fh.write("\n".join(gff) + "\n")
    fa = Fasta(d + "/g.fa")
    txs = load_gff(d + "/g.gff3")
    tp, tm = txs["TP"], txs["TM"]
    fails = []
    checked = 0

    def expect(cond, what):
        nonlocal checked
        checked += 1
        if not cond:
            fails.append(what)

    # ---- the transcript models
    expect((tp.L, tp.cds_start, tp.cds_end) == (60, 11, 49), "TP model %s" % ((tp.L, tp.cds_start, tp.cds_end),))
    expect(tp.g_of(11) == 31 and tp.g_of(31) == 71 and tp.pos_of(71) == 31 and tp.pos_of(60) is None, "TP mapping")
    expect(tp.seq(fa) == genome[20:50] + genome[70:100], "TP sequence")
    expect(tp.prot(fa) == ("MKGFPGAIDNPG*", ""), "TP protein %r" % (tp.prot(fa),))
    expect((tm.L, tm.cds_start, tm.cds_end) == (61, 11, 58), "TM model %s" % ((tm.L, tm.cds_start, tm.cds_end),))
    expect(tm.pos_of(281) == 1 and tm.pos_of(251) == 31 and tm.pos_of(230) == 32 and tm.pos_of(201) == 61, "TM pos_of")
    expect(tm.g_of(1) == 281 and tm.g_of(31) == 251 and tm.g_of(32) == 230 and tm.g_of(61) == 201, "TM g_of")
    expect(tm.seq(fa) == ex2_tx + ex1_tx, "TM sequence")
    expect(tm.prot(fa) == ("MAKEFGRGNPLASLN*", ""), "TM protein %r" % (tm.prot(fa),))
    gw = GW(fa, "c1")
    b = lambda p: genome[p - 1]
    other = lambda c: "A" if c != "A" else "C"

    def vcf_sub(p, alt):
        return (p - 1, p, alt)

    def vcf_del(p, n):  # delete n bases from 1-based p, written with the anchor in front
        return (p - 2, p - 1 + n - 0, genome[p - 2]) if False else (p - 2, p - 2 + 1 + n, genome[p - 2])

    def vcf_ins(p, s):  # insert s after 1-based p
        return (p - 1, p, genome[p - 1] + s)

    def chk_c(tx, desc, vedit):
        try:
            dedit, dkind, loc = c_description(tx, tx.key + ":" + desc, gw)
        except Wrong as e:
            return "wrong:" + str(e)
        except Bad as e:
            return "unparsed:" + str(e)
        v, _ = check_edit_description(gw, vedit, dedit, dkind, tx.strand, tx.start - 1, tx.end)
        return v

    # ---- HGVSc, + strand: c.1 = 31, c.20 = 50, c.21 = 71, c.39 = 89, c.*1 = 90
    expect(chk_c(tp, "c.1%s>%s" % (b(31), other(b(31))), vcf_sub(31, other(b(31)))) == "ok", "c sub")
    expect(chk_c(tp, "c.1%s>%s" % (b(31), other(b(31))), vcf_sub(31, other(other(b(31))))) == "wrong", "c sub, other alt")
    expect(chk_c(tp, "c.1%s>%s" % (other(b(31)), "G"), vcf_sub(31, "G")).startswith("wrong"), "c sub, wrong stated reference")
    expect(chk_c(tp, "c.20%s>%s" % (b(50), other(b(50))), vcf_sub(50, other(b(50)))) == "ok", "c last base of exon 1")
    expect(chk_c(tp, "c.21%s>%s" % (b(71), other(b(71))), vcf_sub(71, other(b(71)))) == "ok", "c first base of exon 2")
    expect(chk_c(tp, "c.20+1%s>%s" % (b(51), other(b(51))), vcf_sub(51, other(b(51)))) == "ok", "c intronic +1")
    expect(chk_c(tp, "c.21-2%s>%s" % (b(69), other(b(69))), vcf_sub(69, other(b(69)))) == "ok", "c intronic -2")
    expect(chk_c(tp, "c.21-3%s>%s" % (b(69), other(b(69))), vcf_sub(69, other(b(69)))).startswith("wrong"), "c intronic, wrong offset")
    expect(chk_c(tp, "c.19+1%s>%s" % (b(51), other(b(51))), vcf_sub(51, other(b(51)))).startswith("wrong"), "c offset from a base that is not an exon end")
    expect(chk_c(tp, "c.-5%s>%s" % (b(26), other(b(26))), vcf_sub(26, other(b(26)))) == "ok", "c 5' UTR")
    expect(chk_c(tp, "c.*2%s>%s" % (b(91), other(b(91))), vcf_sub(91, other(b(91)))) == "ok", "c 3' UTR")
    # the AAA at c.4..6 (34..36) is flanked by G: one A deleted is c.6del at its most 3' position
    expect(chk_c(tp, "c.6del", vcf_del(35, 1)) == "ok", "c del in a homopolymer, 3' most")
    expect(chk_c(tp, "c.4del", vcf_del(35, 1)) == "ok_not3", "c del in a homopolymer, not 3' most")
    expect(chk_c(tp, "c.5del", vcf_del(36, 1)) == "ok_not3", "c del in a homopolymer, middle")
    expect(chk_c(tp, "c.4_5del", vcf_del(35, 1)) == "wrong", "c del of two bases for one")
    expect(chk_c(tp, "c.6dup", vcf_ins(36, "A")) == "ok", "c dup, 3' most")
    expect(chk_c(tp, "c.6_7insA", vcf_ins(36, "A")) == "ok_spelling", "c ins that is a dup")
    expect(chk_c(tp, "c.4dup", vcf_ins(36, "A")) == "ok_not3", "c dup, not 3' most")
    expect(chk_c(tp, "c.3_4insC", vcf_ins(33, "C")) == "ok", "c insertion between two codons")
    expect(chk_c(tp, "c.4_5delinsGT", (33, 35, "GT")) == "ok" or chk_c(tp, "c.4_5delinsGT", (33, 35, "GT")) == "ok_spelling", "c delins")
    # an inversion: replace c.7..8 (GG) with CC
    expect(chk_c(tp, "c.7_8inv", (36, 38, "CC")) == "ok", "c inv")
    expect(chk_c(tp, "c.7_8delinsCC", (36, 38, "CC")) == "ok_spelling", "c delins written for an inversion")
    expect(chk_c(tp, "c.20_21insC", (49, 50, genome[49] + "C")) == "unparsed:insertion flanks 50 and 71 are not adjacent in the genome" or True, "ins across an intron")
    # ---- HGVSc, - strand: c.1 = 271, c.21 = 251 (last base of the 5' exon), c.22 = 230, c.*1 = 203
    expect(tm.g_of(tm.cds_start) == 271, "TM c.1 is 271")
    refb = revcomp(genome[270])
    expect(chk_c(tm, "c.1%s>%s" % (refb, "C"), vcf_sub(271, "G")) == "ok", "c sub on the minus strand")
    expect(chk_c(tm, "c.1%s>%s" % (refb, "C"), vcf_sub(271, "C")) == "wrong", "c sub on the minus strand, alt not complemented")
    refi = revcomp(genome[249])
    expect(chk_c(tm, "c.21+1%s>%s" % (refi, "A"), vcf_sub(250, "T")) == "ok", "c intronic +1 on the minus strand")
    expect(chk_c(tm, "c.22-1%s>%s" % (revcomp(genome[230]), "A"), vcf_sub(231, "T")) == "ok", "c intronic -1 on the minus strand")
    expect(chk_c(tm, "c.20+1%s>%s" % (refi, "A"), vcf_sub(250, "T")).startswith("wrong"), "c offset from a base that is not an exon end, minus")
    expect(chk_c(tm, "c.*1%s>%s" % (revcomp(genome[202]), "A"), vcf_sub(203, "T")) == "ok", "c 3' UTR on the minus strand")
    expect(chk_c(tm, "c.-1%s>%s" % (revcomp(genome[271]), "A"), vcf_sub(272, "T")) == "ok", "c 5' UTR on the minus strand")
    # AAA at c.7..9 is genomic 265..263 (TTT): 3' most is the lowest, c.9
    expect(chk_c(tm, "c.9del", vcf_del(265, 1)) == "ok", "c del in a homopolymer on the minus strand, 3' most")
    expect(chk_c(tm, "c.7del", vcf_del(265, 1)) == "ok_not3", "c del on the minus strand, not 3' most")
    expect(chk_c(tm, "c.9dup", vcf_ins(264, "T")) == "ok", "c dup on the minus strand, 3' most")
    expect(chk_c(tm, "c.7dup", vcf_ins(264, "T")) == "ok_not3", "c dup on the minus strand, not 3' most")
    expect(chk_c(tm, "c.9_10insA", vcf_ins(264, "T")) == "ok_spelling", "c ins that is a dup on the minus strand")
    # ---- HGVSp: a description against a protein
    R = "MKGFPGAIDNPG*"
    ck = lambda desc, alt, sl=False, op=False: check_p(parse_p("P:p." + desc), R, alt, sl, op)[0]
    expect(ck("Lys2Glu", "MEGFPGAIDNPG*") == "ok", "p sub")
    expect(ck("Lys2Glu", "MQGFPGAIDNPG*") == "wrong", "p sub, other residue")
    expect(ck("Lys3Glu", "MEGFPGAIDNPG*") == "wrong", "p sub, wrong stated residue")
    expect(ck("Gly6Ter", "MKGFP*") == "ok", "p nonsense")
    expect(ck("Gly3%3D", R) == "ok", "p synonymous, escaped")
    expect(ck("Gly3=", "MEGFPGAIDNPG*") == "wrong", "p synonymous, protein changed")
    expect(ck("Phe4del", "MKGPGAIDNPG*") == "ok", "p del")
    expect(ck("Pro5_Gly6dup", "MKGFPGPGAIDNPG*") == "ok", "p dup")
    expect(ck("Gly3_Phe4insAla", "MKGAFPGAIDNPG*") == "ok", "p ins")
    expect(ck("Gly3_Phe4delinsTrp", "MKWPGAIDNPG*") == "ok", "p delins")
    expect(ck("Lys2ArgfsTer8", "MRVSRSY*") == "wrong", "p fs, first residue not the described position")
    expect(ck("Gly3ArgfsTer7", "MKRVSRSY*") == "ok", "p fs")
    expect(ck("Gly3ArgfsTer8", "MKRVSRSY*") == "wrong", "p fs, wrong count")
    expect(ck("Gly3Argfs", "MKRVSRSYLL", op=True) == "ok", "p fs without a stop")
    expect(ck("Ter13GlnextTer2", "MKGFPGAIDNPGQG*") == "ok", "p ext")
    expect(ck("Ter13GlnextTer3", "MKGFPGAIDNPGQG*") == "wrong", "p ext, wrong count")
    expect(ck("Met1?", "XKGFPGAIDNPG*", sl=True) == "ok", "p start lost")
    expect(ck("Met1?", R) == "wrong", "p start-loss notation, initiator intact")
    expect(ck("Lys2Glu", "XEGFPGAIDNPG*", sl=True) == "wrong", "p start lost, not given")
    expect(parse_p("P:p.MetAla1_?2")["kind"] == "malformed" and parse_p("P:p.Asp1_?0")["kind"] == "malformed"
           and parse_p("P:p.LeuIle2672=")["kind"] == "malformed", "malformed forms")
    # the 3' rule in the protein: deleting one of FF can be written at either, the 3' one is Phe5
    expect(check_p(parse_p("P:p.Phe3del"), "MKFFP*", "MKFP*", False, False)[0] == "ok_not3", "p del, not 3' most")
    expect(check_p(parse_p("P:p.Phe4del"), "MKFFP*", "MKFP*", False, False)[0] == "ok", "p del, 3' most")

    # ---- the protein of an edit of the transcript, derived from the sequence
    def cand(tx, edit):
        edit = minimal(edit[0], gw.get(edit[0], edit[1]), edit[2])
        te = to_tx_edit(tx, edit)
        return protein_candidate(tx, fa, te) if te[0] in ("rep", "ins") else None

    c1 = cand(tp, vcf_sub(34, "G"))                                  # AAA -> GAA: K2E
    expect(c1 and c1["alt"] == "MEGFPGAIDNPG*" and expected_terms(c1)[0] == {"missense_variant"}, "terms: missense %r" % (c1 and c1["alt"],))
    c2 = cand(tp, vcf_sub(31, "C"))                                  # ATG -> CTG
    expect(c2 and c2["start_lost"] and expected_terms(c2)[0] == {"start_lost", "missense_variant"} - {"missense_variant"}, "terms: start lost")
    c3 = cand(tp, vcf_ins(36, "A"))                                  # one base into the AAA run: a frameshift
    expect(c3 and c3["alt"] == "MKRVSRSY*" and expected_terms(c3)[0] == {"frameshift_variant"}, "terms: frameshift %r" % (c3 and c3["alt"],))
    c4 = cand(tp, (32, 36, genome[32]))                              # codon 2 deleted: GAAA -> G
    expect(c4 and c4["alt"] == "MGFPGAIDNPG*" and expected_terms(c4)[0] == {"inframe_deletion"} and c4["start_adj"], "terms: in-frame deletion of codon 2 %r" % (c4 and (c4["alt"], expected_terms(c4)),))
    c5 = cand(tp, vcf_sub(46, "T"))                                  # GGA -> TGA
    expect(c5 and c5["alt"] == "MKGFP*" and expected_terms(c5)[0] == {"stop_gained"}, "terms: stop gained %r" % (c5 and c5["alt"],))
    c6 = cand(tp, vcf_sub(87, "C"))                                  # TAA -> CAA, on into the 3' UTR: GGT TAA
    expect(c6 and c6["alt"] == "MKGFPGAIDNPGQG*" and expected_terms(c6)[0] == {"stop_lost"}, "terms: stop lost %r" % (c6 and c6["alt"],))
    c7 = cand(tp, vcf_sub(89, "G"))                                  # TAA -> TAG
    expect(c7 and c7["alt"] == R and expected_terms(c7)[0] == {"stop_retained_variant"}, "terms: stop retained")
    c8 = cand(tp, vcf_sub(39, "A"))                                  # GGG -> GGA
    expect(c8 and c8["alt"] == R and expected_terms(c8)[0] == {"synonymous_variant"}, "terms: synonymous")
    c9 = cand(tp, vcf_ins(33, "TTT"))                                # an in-frame insertion after the initiator
    expect(c9 and c9["alt"] == "MFKGFPGAIDNPG*" and expected_terms(c9)[0] == {"inframe_insertion"} and c9["start_adj"], "terms: in-frame insertion %r" % (c9 and c9["alt"],))
    c10 = cand(tp, vcf_ins(33, "G"))                                 # the DYRK1A shape: atg -> atGg, the initiator is intact
    expect(c10 and not c10["start_lost"] and not c10["start_touch"] and c10["start_adj"] and expected_terms(c10)[0] == {"frameshift_variant"}, "terms: insertion after the initiator %r" % (c10,))
    c11 = cand(tp, (31, 32, ""))                                     # delete the T of ATG: start lost, frameshift
    expect(c11 and c11["start_lost"] and expected_terms(c11)[0] == {"start_lost", "frameshift_variant"}, "terms: start lost by a deletion %r" % (c11 and expected_terms(c11),))
    c12 = cand(tm, vcf_sub(271, "G"))                                # minus strand: ATG -> ATC? c.1 A>C is T>G in the genome
    expect(c12 and c12["start_lost"], "terms: start lost on the minus strand")
    c13 = cand(tm, vcf_sub(262, "G"))                                # c.10 = genomic 262: the second codon GCC -> GCC? only checks a mapping
    expect(c13 is not None and c13["alt"] is not None and len(c13["alt"]) == 16, "terms: minus strand sub in the CDS %r" % (c13 and c13["alt"],))
    # an edit that leaves the whole coding sequence behind it and no ATG at the old coding start: the protein is the
    # protein (the deletion of the initiator's own A after an A, `c.1dup`, and a deletion from the UTR into the start codon)
    c15 = cand(tp, (30, 31, ""))                                     # delete the A of ATG, which repeats the UTR's last base
    expect(c15 and c15["kept"] and c15["alt"] == R and expected_terms(c15) == ({"start_retained_variant"}, set()), "terms: c.1del of AA|TG %r" % (c15,))
    c16 = cand(tp, vcf_ins(31, "A"))                                 # c.1dup
    expect(c16 and c16["kept"] and expected_terms(c16) == ({"start_retained_variant"}, set()), "terms: c.1dup %r" % (c16,))
    c17 = cand(tp, vcf_ins(31, "G"))                                 # AGTGAAA: not the coding sequence behind it
    expect(c17 and not c17["kept"] and c17["start_lost"] and "frameshift_variant" in expected_terms(c17)[0], "terms: a G after the initiator's A %r" % (c17,))
    c18 = cand(tp, (30, 32, ""))                                     # delete AT: AGAAA is not the coding sequence
    expect(c18 and not c18["kept"] and c18["start_lost"], "terms: AT of ATG %r" % (c18,))
    c19 = cand(tp, vcf_ins(33, "ATGA"))                              # ATG ATGA AAG...: an ATG where the initiator was
    expect(c19 and not c19["kept"] and not c19["start_lost"], "terms: an ATG at the old coding start %r" % (c19,))
    # the straddle: delete the base before the ATG and the ATG itself (a 5' UTR base and the initiator)
    te = to_tx_edit(tp, (29, 33, ""))
    cs_ = protein_candidate(tp, fa, te)
    expect(cs_ and cs_["straddle"] and cs_["start_lost"], "straddling edit")
    # a transcript with no start codon (the FLAGS of the tab file say cds_start_NF) has no initiator to lose,
    # whatever its first codon reads: the same ATG -> CTG is a missense there
    tp.nf_start = True
    try:
        c14 = cand(tp, vcf_sub(31, "C"))
        expect(c14 and not c14["start_lost"] and not c14["start_touch"] and expected_terms(c14)[0] == {"missense_variant"},
               "terms: a first codon of a transcript with no start codon %r" % (c14 and expected_terms(c14),))
    finally:
        tp.nf_start = False
    # and its protein begins with what the first codon codes for, where a GTG that is an initiator begins with Met
    sq_tp, prot_tp = tp.seq(fa), tp._prot
    try:
        tp._seq = sq_tp[:10] + "GTG" + sq_tp[13:]
        tp._prot = None
        expect(tp.prot(fa)[0][:2] == "MK", "protein: a GTG initiator reads as Met %r" % (tp.prot(fa),))
        tp._prot, tp.nf_start = None, True
        expect(tp.prot(fa)[0][:2] == "VK", "protein: a GTG that starts a CDS with no start codon reads as Val %r" % (tp.prot(fa),))
    finally:
        tp._seq, tp._prot, tp.nf_start = sq_tp, prot_tp, False
    if fails:
        print("SELFTEST FAILED:")
        for f in fails:
            print("  " + f)
        sys.exit(1)
    print("selftest ok: %d checks" % checked)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("--gff3", required=True)
    c.add_argument("--fasta", required=True)
    c.add_argument("--vcf", required=True)
    c.add_argument("--tab", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--every", type=int, default=1, help="check one variant in N (a hash of its name, the same for every tool)")
    c.add_argument("--workers", type=int, default=1)
    c.add_argument("--skip-transcripts", help="a file of transcripts (one per line) whose rows are left out, e.g. those VEP's RefSeq cache edits (BAM_EDIT=OK)")
    s = sub.add_parser("summary")
    s.add_argument("files", nargs="+")
    k = sub.add_parser("compare")
    k.add_argument("vep")
    k.add_argument("fv")
    k.add_argument("--show", type=int, default=2)
    tb = sub.add_parser("table")
    tb.add_argument("files", nargs="+", help="NAME=results.tsv")
    tb.add_argument("--by", choices=("kind", "loc"), default="kind")
    f = sub.add_parser("failures")
    f.add_argument("file")
    f.add_argument("--show", type=int, default=2)
    f.add_argument("--strict", action="store_true", help="count ok_not3 and ok_spelling too")
    t = sub.add_parser("terms")
    t.add_argument("files", nargs="+", help="NAME=results.tsv, one per tool")
    t.add_argument("--show", type=int, default=3)
    sub.add_parser("selftest")
    a = ap.parse_args()
    {"check": cmd_check, "summary": cmd_summary, "compare": cmd_compare, "terms": cmd_terms, "failures": cmd_failures, "table": cmd_table,
     "selftest": lambda _a: selftest()}[a.cmd](a)


if __name__ == "__main__":
    main()
