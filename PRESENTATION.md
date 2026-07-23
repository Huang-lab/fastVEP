# echtvar: Efficient Compact Compressed Haplotype and Variant Archive

---

## Slide 1: The Problem

**Goal:** Efficiently store and query annotations for **~800 million variants** with hundreds of possible INFO fields.

### The scale

- **~800M variants** in the dataset
- **~11 commonly queried annotation fields per variant** (AC, AN, AF, nhomalt, FILTER, etc.)
- **~8.8B annotation values** across those fields
- Must support fast lookups:
  - **Query:** `(chromosome, position, ref, alt)`
  - **Result:** all annotations for that variant
- Must compress well: raw archives are **tens to hundreds of GB**

### Naive approach: Row-oriented storage

Store each variant as a row, with each annotation as a column:

- **800M rows × ~11 fields = ~8.8B values**
- If values average just **4 bytes** → **~35 GB** of raw numeric data
- Add variant identity, offsets, null/missing-value markers, strings, and metadata → **potentially 50–100+ GB uncompressed**
- General-purpose compression helps, but random point lookups can be expensive
- Indexing every field adds additional storage and complexity

### echtvar's insight

**Pack variant identity into a single sortable 32-bit integer**, allowing variants to be stored in sorted, compressed arrays.

This enables:

- **Binary search** to locate a variant
- **Compact columnar arrays** for annotations
- **Efficient compression** by exploiting sorted genomic positions and repeated values
- **No large per-field database index required**

> **Key idea:** Turn variant lookup into a binary search over a compact, sorted representation rather than scanning a massive table.

## Slide 2: var32 — The Core Encoding

**A 32-bit bitfield that encodes a short variant (≤4 bases total) into one sortable key:**

```
u32 layout (LSB → MSB):
┌─────────────┬──────┬──────┬────────────┐
│  enc (8b)   │alen  │ rlen │  position  │
│   0-7       │ 8-9  │10-11 │  12-31     │
└─────────────┴──────┴──────┴────────────┘
```

**Each field:**
- `position` (20 bits): genomic position within a 1MB chunk (0..1,048,575)
- `rlen` (2 bits): ref allele length (1-3, or 3=too-long flag)
- `alen` (2 bits): alt allele length (1-3, or 3=too-long flag)
- `enc` (8 bits): ref+alt bases 2-bit-packed (A=0, C=1, G=2, T=3)

**Example:** Position 423432, REF="A", ALT="ACA"
```
rlen=1, alen=3
enc encoding (left-to-right: REF then ALT):
  A (0b00) → 0
  A (0b00) → 0
  C (0b01) → 1
  A (0b00) → 4   (shift and accumulate: 0*4+0, 0*4+0, 0*4+1, 1*4+0)
Final: u32 = 1,734,379,268
```

**Why it works:**
- Sorts naturally: highest bits (position) dominate, then rlen/alen, then enc
- Binary-searchable: no custom comparator needed
- Dense: 4 bases + metadata in 32 bits
- Supports ref/alt boundary detection: rlen/alen stored separately, not inferred from enc

---

## Slide 3: kmer16 — Handling Long Variants

**For variants where ref.len() + alt.len() > 4 (too long for var32):**

Uses the same 2-bit-per-base encoding but spills into a `Vec<u32>`:

```rust
pub struct LongVariant {
    position: u32,
    idx: u32,
    sequence: Vec<u32>,   // kmer16 encoded
}
```

**Encoding:** Same 2-bit packing but packed into multiple u32s (16 bases per word):
```
sequence[0] = ref_len (full u32)
sequence[1] = alt_len (full u32)
sequence[2..n] = packed bases (ref then alt, 16 bases per u32)
```

**Example:** REF="AAAAC", ALT="A" (6 bases total, needs 3 words)
```
sequence[0] = 5         (ref length)
sequence[1] = 1         (alt length)
sequence[2] = 4         (A=0,A=0,A=0,A=0,C=1,A=0 packed into one u32)
```

**Lookup:** Same binary search as var32, but compares `Vec<u32>` lexicographically
- Sorted by (position, sequence)
- Stored separately in `too-long-for-var32.enc` within each chunk
- Same parallel value arrays indexed by `idx` field

---

## Slide 4: Quantization — Storing Floats as Integers

**Problem:** Stream-vbyte compression (the core of echtvar's pipeline) only handles integers efficiently.

**Solution:** Store floats as scaled integers; multiply on write, divide on read.

```
Float value: 0.0000015 (e.g., allele frequency)
Multiplier:  2,000,000 (from config.json)
On write:    0.0000015 × 2,000,000 = 3 (rounded to u32)
On disk:     3 (compresses as small integer)
On read:     3 / 2,000,000 = 0.0000015 (recovered)
```

**From gnomAD config:**
```json
{"field": "AF", "alias": "gnomad_af", "multiplier": 2000000}
```

**Why it works:**
- Keeps everything in the u32 pipeline (same stream-vbyte encoding as integer fields)
- Small multipliers (1–10 million) trade quantization error for compression efficiency
- Quantization error ≤ 1/multiplier (e.g., ±0.0000005 for AF with multiplier 2M)

**Combined with zigzag (for signed values):**
```rust
fn get_float_value(idx: usize) -> f32 {
    let v: u32 = raw_value[idx];
    if v == u32::MAX {
        missing_value
    } else if field.zigzag {
        (zigzag::decode(v) as f32) / field.multiplier as f32
    } else {
        v as f32 / field.multiplier as f32
    }
}
```

---

## Slide 5: Stream-vbyte + Delta Encoding — Compression

**Layer 1: Delta encoding** (because var32s are sorted, deltas are small)
```
Original:  [100, 102, 105, 110, ...]
Delta:     [100, 2, 3, 5, ...]  (each value is diff from previous)
Benefit:   Small numbers compress better
```

**Layer 2: Stream-vbyte** (variable-length integer codec)
- Encodes a control stream (2 bits per value: how many bytes this value takes)
- Data stream (packed variable-length encoded bytes)
- Uses SIMD (SSSE3 shuffles on x86_64) to decode ~100M integers/sec

```
Control bytes: | 2 | 1 | 3 | 2 | ...  (says: next int takes 2 bytes, then 1, then 3, then 2)
Data bytes:    | byte0 byte1 | byte0 | byte0 byte1 byte2 | byte0 byte1 | ...
```

**Compression ratio:** ~4:1 on real gnomAD data (var32 column compresses to ~1/4 original size)

**Decoding cost:** ~once per chunk, not per variant
- On `set_position` (chunk transition), entire chunk's var32s decoded into memory at once
- Cumulative sum applied to reconstruct from deltas
- Then reused for all lookups within that chunk

---

## Slide 6: Chunking — Spatial Organization

**File layout:** Chunks are the horizontal partitioning unit

```
echtvar/
├── config.json                    (field metadata)
├── strings/
│   └── gnomad_filter.txt          (categorical lookup table)
└── {chrom}/
    ├── 0/                         (chunk 0: positions 0..1,048,575)
    │   ├── var32.bin              (delta + stream-vbyte encoded var32s)
    │   ├── gnomad_ac.bin          (stream-vbyte encoded values)
    │   ├── gnomad_an.bin
    │   ├── gnomad_af.bin
    │   └── too-long-for-var32.enc (bincode-serialized LongVariants)
    ├── 1/                         (chunk 1: positions 1,048,576..2,097,151)
    │   ├── var32.bin
    │   ├── gnomad_ac.bin
    │   └── ...
    └── ...
```

**Why 1<<20 (1,048,576)?**
- Fits in 20 bits (the `position` field of var32)
- Chunk index = `pos >> 20` (fast division)
- Within-chunk offset = `pos & 0xFFFFF` (fast modulo)
- ~2,959 chunks for whole genome (GRCh38)

**Lookup:** To find a variant at position P:
1. Open `echtvar/{chrom}/{P >> 20}/`
2. Decode all fields for that chunk (once, cached)
3. Binary search sorted var32s
4. Index into value arrays

**Benefit:** Chunk boundary is the only place where expensive decode happens; amortized across thousands of variants.

---

## Slide 7: The Lookup Pipeline — End-to-End

**Query:** Given a variant at runtime (chrom, pos, ref, alt), retrieve all annotations

**Step 1: Encode query variant**
```rust
let query_var32 = var32::encode(pos, ref_bytes, alt_bytes);
// Produces same 32-bit key as if it were stored
```

**Step 2: Determine chunk**
```rust
let chunk_idx = pos >> 20;
if chunk_idx != last_chunk {
    // Load new chunk: decode all field columns
    echtvar.set_position(rid, chrom, pos);
}
```

**Step 3: Binary search**
```rust
match echtvar.var32s.binary_search(&query_var32) {
    Ok(idx) => {
        // Found: index into parallel value arrays
        for field in fields {
            let value = echtvar.values[field.values_i][idx];
            // Decode: handle missing sentinel, zigzag, multiplier
        }
    }
    Err(_) => {
        // Not found: use missing values
    }
}
```

**Complexity:**
- Per-chunk: O(chunk_size) decode (amortized across variants in chunk)
- Per-variant: O(log N) binary search + O(1) value indexing

**Performance:** ~2.0M variants/sec on sorted VCF, single-threaded

---

## Slide 8: Missing Sentinel — u32::MAX

**Problem:** How to represent "this variant has no value for this field"?

**Solution:** Reserve the single largest u32 value as a sentinel.

```rust
const MISSING: u32 = u32::MAX;  // 4,294,967,295
```

**On read:**
```rust
fn get_int_value(idx: usize) -> i32 {
    let v: u32 = raw_values[idx];
    if v == u32::MAX {
        field.missing_value  // Return configured missing value (e.g., -1)
    } else {
        v as i32  // Normal value
    }
}
```

**Cost:** Every value lookup branches on this check; can't use full u32 range

**Alternative (unused here): Arrow's approach**
- Separate **validity bitmap** (1 bit per value, 0=null)
- Data buffer holds garbage at null slots (not checked)
- Allows using full u32 range for real data

**Why echtvar chose sentinel:**
- Simpler, no separate bitmap
- Stream-vbyte/delta pipeline works on flat u32 arrays
- Branch is cheap and sparse (most variants have most fields)

---

## Slide 9: Zigzag Encoding — Handling Negative Values

**Problem:** Stream-vbyte compresses small unsigned integers well. But some fields (e.g., CADD raw scores) can be negative.

**Solution:** Map signed integers to unsigned so small-magnitude negatives stay small.

```
Standard cast: -1_i32 as u32 = 4,294,967,295 (huge, collides with MISSING sentinel!)
Zigzag:        -1_i32 → 1_u32 (small, compresses well)

Mapping:  ... -2 → 3, -1 → 1, 0 → 0, 1 → 2, 2 → 4, ...
```

**Zigzag encode/decode:**
```rust
fn zigzag_encode(i: i32) -> u32 {
    ((i << 1) ^ (i >> 31)) as u32
}
fn zigzag_decode(e: u32) -> i32 {
    ((e >> 1) as i32) ^ -((e & 1) as i32)
}
```

**Per-field opt-in:** Config file specifies `zigzag: true` for signed fields
```json
[
  {"field": "CADD_raw", "alias": "cadd_score", "zigzag": true}
]
```

**On read (combined with quantization):**
```rust
if field.zigzag {
    (zigzag::decode(v) as f32) / field.multiplier as f32
} else {
    (v as f32) / field.multiplier as f32
}
```

---

## Slide 10: Categorical String Tables — Efficient Encoding

**Problem:** Categorical fields (like FILTER="PASS/LOW_CONFIDENCE") are repetitive; storing strings inline wastes space.

**Solution:** Assign integers to unique strings; store integers, keep a lookup table.

```
String values:      "PASS", "PASS", "LOW_CONF", "PASS", ...
Indexed:            0, 0, 1, 0, ...
Lookup table:       strings[0]="PASS", strings[1]="LOW_CONF"
Compressed size:    ~1/10th of storing strings inline
```

**Storage:**
```
echtvar/
├── config.json
└── strings/
    └── gnomad_filter.txt    (one string per line)
        PASS
        LOW_CONFIDENCE
        MISSING
```

**Field metadata (config.json):**
```json
{
  "field": "FILTER",
  "alias": "gnomad_filter",
  "ftype": "Categorical",
  "missing_string": "MISSING"
}
```

**On write:** Map string → index, store index as u32 (same pipeline as integers)

**On read:**
```rust
fn get_categorical_value(idx: usize) -> String {
    let i: u32 = raw_values[idx];
    if i == u32::MAX {
        field.missing_string.clone()
    } else {
        field_strings[i as usize].clone()
    }
}
```

**Cardinality limit:** Only expose strings to filter expressions if <256 unique values (otherwise evaluation becomes slow)

---

## Slide 11: Sorted-VCF Workload — Why the Architecture Wins

**Typical annotation workflow:** Stream a sorted VCF, annotate each variant, output results
```
Input VCF: sorted by (chrom, position)
echtvar:   organized by (chrom, 1MB chunks)
→ Natural alignment: variants arrive in chunk order
```

**Why this matters:**

**Sequential access pattern:**
```
VCF stream:   chr1:100 → chr1:500 → chr1:2000 → chr1:50000 → ...
Chunks hit:   chunk 0  → chunk 0  → chunk 0   → chunk 0    → ...
Decode cost:  paid once, amortized across ~1M variants per chunk
```

**Binary search + indexing (per-variant):**
- O(log N) comparisons on already-resident, cache-warm buffer
- Small constant work, hits CPU cache well

**Result:** ~2.0M variants/sec single-threaded

**Contrast: random access pattern**
```
Query 1: chr5:12345
Query 2: chr3:98765
Query 3: chr22:1000000
Chunks:  each query hits different chunk
Decode:  paid 3 times, rarely amortized
Result:  ~20x slower per query
```

**Lesson:** echtvar is optimized for its common use case (annotating a sorted VCF); random ad-hoc queries pay full per-chunk decode cost.

---

## Slide 12: The Parquet Alternative — What Changed, What Didn't

**Proposal:** Swap zip container for Apache Parquet columnar format. Keep variant encoding identical.

**What stays the same:**
- var32 encoding (same 32-bit key, same 20-bit position field)
- Long-variant path (kmer16 encoding)
- Quantization + multiplier (floats as scaled integers)
- Zigzag for signed values
- Missing sentinel (u32::MAX)
- Chunking strategy (one Parquet row group per 1MB chunk)

**What changes:**
```
zip layout:                           Parquet layout:
echtvar/chr1/0/var32.bin             One .parquet file with:
echtvar/chr1/0/gnomad_ac.bin     →   ├── row group 0 (chunk 0)
echtvar/chr1/0/gnomad_an.bin         │   ├── column: var32 (sorted)
echtvar/chr1/0/gnomad_af.bin         │   ├── column: gnomad_ac
... ~38,500 zip entries total         │   ├── column: gnomad_an
config.json + strings/* separate      │   └── ...
                                      ├── row group 1 (chunk 1)
                                      │   └── ...
                                      └── single footer (schema + all column metadata)
```

**Key findings:**

1. **Size:** ~1.12x vs zip (with delta-encoding fix on var32 column)
2. **Sequential throughput:** ~2.56M variants/sec vs ~2.0M (1.28x faster)
   - Per-variant search: 6.3M lookups/sec (3.2x faster than zip's combined rate)
   - But decode cost was also measured per-chunk
3. **Correctness:** All ~13M lookups from gnomAD matched bit-for-bit

**The real win:** Standard columnar format readable by DuckDB/Polars/pandas without echtvar code — enables ad-hoc querying of encoded population data

**The cost:** Parquet + Arrow are heavier dependencies than zip + stream-vbyte

**Caveats:**
- Measured as standalone microbenchmark (binary_search + column index), not full htslib pipeline
- Requires parsing footer once per session (naive re-parse per chunk = 2x slower than zip)
- Random access still pays full per-chunk decode (same as zip)

---

## Summary: Design Principles

1. **Pack variant identity efficiently** (var32) to enable binary search
2. **Keep everything integer-shaped** (quantization, zigzag) for compression
3. **Amortize expensive operations** (per-chunk decode, cumulative sum) across many queries
4. **Optimize for sorted workloads** (sequential access, chunk locality)
5. **Separate concerns:** encoding (var32/kmer16) from compression (stream-vbyte) from organization (chunks)

**Result:** ~13M variants + ~11 fields compressed to ~1.1x the VCF size, searchable in 0.5µs per variant on sorted VCF.

---
---

# fastVEP: Effect Prediction

**The problem:** given (pos, ref, alt) and a set of overlapping transcripts, compute what happens biologically.
→ A geometry + biology problem. The hard part is coordinate transforms, reading frames, and splice/protein arithmetic.

```
variant → [map onto transcript] → [translate codons] → consequence (derived)
```

This requires an entirely dedicated core: transcript models, coordinate maps, codon tables, and a severity ontology.

---

## Slide 13: The Consequence Ontology — 49 Ranked Terms, 4 Impact Tiers

**Every predicted effect is one of 49 Sequence-Ontology terms**, each with a numeric rank (1 = most severe) and a coarse `Impact`:

```
Impact:      High  >  Moderate  >  Low  >  Modifier
Examples:    transcript_ablation (rank 1, High)
             missense_variant    (rank 12, Moderate)
             synonymous_variant  (rank 20, Low)
             intron_variant      (rank 38, Modifier)
```

**Why it matters:** a variant can overlap several transcripts (or isoforms), each producing a different consequence.

```rust
fn most_severe(consequences: &[Consequence]) -> Consequence {
    consequences.iter().min_by_key(|c| c.rank()).copied().unwrap()
}
```

This ontology is the backbone of every downstream report.

`crates/fastvep-core/src/consequence.rs`

---

## Slide 14: Coordinate Transform Pipeline — Genomic → cDNA → CDS → Protein

**fastVEP walks through four coordinate systems** to know what a base change *means*:

```
genomic position  →  cDNA position  →  CDS position  →  protein (codon) position
  (chr:12345)          (c.152)           (subtract UTR)      (152-1)/3 + 1 = 51
```

**genomic_to_cdna:** walk sorted exons, accumulate spliced length, flip direction on minus strand
**cdna_to_cds:** subtract coding-region start, add reading-frame phase
**cds_to_protein:** pure arithmetic — `(cds_pos - 1) / 3 + 1`

```rust
// closed-form, not a lookup table — recomputed per query
let cds_pos = cdna_pos - coding_start + phase;
let protein_pos = (cds_pos - 1) / 3 + 1;
```

`crates/fastvep-genome/src/transcript.rs:94-127`

## Slide 15: Concurrency Model — Batched Parallel Pipeline

**fastVEP runs a three-phase batched pipeline per VCF chunk:**
```
Phase 1    (sequential):  read/buffer a batch of VCF records
Phase 1.5  (sequential):  preload/prime transcript + SA caches for the batch's positions
                           (avoids redundant per-thread I/O and cache misses)
Phase 2    (parallel):    par_iter_mut() over the batch — transcript lookup,
                           consequence prediction, HGVS generation, all per-thread
```

**Why this shape:** effect prediction is CPU-heavy (coordinate transforms, codon translation, HGVS string building) per variant — genuinely parallelizable work.

**Result:** 46k–86k variants/sec depending on organism; 2.6x–130x speedup over Ensembl VEP.

`crates/fastvep-cli/src/pipeline.rs:600-650`, `README.md:573-600`
