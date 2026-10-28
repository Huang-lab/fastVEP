#!/usr/bin/env bash
# Rebuilds the binary fixtures in this directory. Needs bgzip and tabix (htslib).
#
# The records are synthetic. Their shapes are the ones issue #139 reports: an
# anchored deletion and an anchored insertion that VEP matches by trimming, a
# minimal deletion and insertion, an SNV site with a multi-allelic, a failed and
# an allele-less record on it, and a mitochondrial record.
#
# The index is made the way VEP's own converter makes it
# (ensembl-vep/convert_cache.pl, process_chr_type):
#
#     bgzip all_vars
#     tabix -C -s 1 -b <start column> -e <start column> all_vars.gz
#
# `-C` writes a .csi, which is why a cache converted by VEP carries no .tbi. The
# start column is both -b and -e, so the index covers each record's start only,
# not its span. `tbi/` is the same data indexed without -C, the form older caches
# have. `noindex/` has the data and no index at all. `endcol/` is indexed on the
# end column instead: tabix accepts that, and it answers queries with the wrong
# records, which is why the reader checks which column an index was built on.
# `both/` has the right .csi and a wrong .tbi side by side, to pin which one is
# read when a directory holds both: the `tabix` command takes the .csi.
#
# info.txt follows what convert_cache.pl leaves behind (a leading `chr` in
# variation_cols, `var_type tabix`, no valid_chromosomes line: VEP derives that
# from the sub-directories, CacheDir.pm). chr_synonyms.txt is hand-written in the
# format CacheDir.pm auto-loads from a cache directory; its contents are
# illustrative, not copied from a real cache.
set -euo pipefail
# byte-order sorting, so records that share a start land in the same order on
# every machine (the tests read that order)
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
cols="chr,variation_name,failed,somatic,start,end,allele_string,strand,minor_allele,minor_allele_freq,clin_sig,phenotype_or_disease,clin_sig_allele,pubmed,var_synonyms,AFR,AMR,EAS,EUR,SAS"
start_col=5 # 1-based, counting the leading chr column

records_10() {
  cat <<'REC'
10|rsFIX_DEL_MIN|0|0|10400|10404|AACCC/-|1|.|.|.|0|.|.|.|.|.|.|.|.
10|rsFIX_DEL_ANCH|0|0|10401|10406|ACCCTA/A|1|.|.|.|0|.|.|.|.|.|.|.|.
10|rsFIX_INS_MIN|0|0|10409|10408|-/T|1|.|.|.|0|.|.|.|.|.|.|.|.
10|rsFIX_SNV|0|0|10414|10414|C/A|1|A|0.0020|benign|0|.|.|.|A:0.0100|.|.|.|.
10|rsFIX_SNV_MULTI|0|0|10414|10414|C/A/G|1|.|.|.|0|.|.|.|.|.|.|.|.
10|rsFIX_SNV_FAILED|1|0|10414|10414|C/A|1|.|.|.|0|.|.|.|.|.|.|.|.
10|COSVFIX1|0|1|10414|10414|COSMIC_MUTATION|1|.|.|.|0|.|.|.|.|.|.|.|.
10|rsFIX_INS_ANCH|0|0|10420|10420|G/GA|1|.|.|.|0|.|.|.|.|.|.|.|.
REC
}

records_mt() {
  cat <<'REC'
MT|rsFIX_MT|0|0|100|100|A/G|1|.|.|.|0|.|.|.|.|.|.|.|.
REC
}

info() {
  printf '%s\n' \
    '# CACHE UPDATED 2026-04-08 00:00:00' \
    'species	homo_sapiens' \
    'assembly	GRCh38' \
    'cache_region_size	1000000' \
    "variation_cols	$cols" \
    'var_type	tabix'
}

# sort exactly as tabix needs: by contig, then by the start column
write_chr() { # <dest dir> <chr> <records fn> <csi|tbi|endcol|both|none>
  local dest="$1" chr="$2" fn="$3" mode="$4"
  mkdir -p "$dest/$chr"
  "$fn" | tr '|' '\t' | sort -t$'\t' -k1,1 -k"$start_col","$start_col"n >"$dest/$chr/all_vars"
  bgzip -f "$dest/$chr/all_vars"
  case "$mode" in
    csi) tabix -C -s 1 -b "$start_col" -e "$start_col" "$dest/$chr/all_vars.gz" ;;
    tbi) tabix -s 1 -b "$start_col" -e "$start_col" "$dest/$chr/all_vars.gz" ;;
    # indexed on the end column instead of start: searchable, and wrong
    endcol) tabix -C -s 1 -b "$((start_col + 1))" -e "$((start_col + 1))" "$dest/$chr/all_vars.gz" ;;
    # the right .csi, then a .tbi on the end column beside it
    both)
      tabix -C -s 1 -b "$start_col" -e "$start_col" "$dest/$chr/all_vars.gz"
      tabix -s 1 -b "$((start_col + 1))" -e "$((start_col + 1))" "$dest/$chr/all_vars.gz"
      ;;
    none) ;;
    *) echo "unknown index mode $mode" >&2; exit 1 ;;
  esac
}

rm -rf "$here/csi" "$here/tbi" "$here/noindex" "$here/endcol" "$here/both"

for d in csi tbi noindex endcol both; do mkdir -p "$here/$d"; info >"$here/$d/info.txt"; done
printf '%s\n' '10	chr10	NC_000010.11' 'MT	chrM	NC_012920.1' >"$here/csi/chr_synonyms.txt"

write_chr "$here/csi" 10 records_10 csi
write_chr "$here/csi" MT records_mt csi
write_chr "$here/tbi" 10 records_10 tbi
write_chr "$here/tbi" MT records_mt tbi
write_chr "$here/noindex" 10 records_10 none
write_chr "$here/endcol" 10 records_10 endcol
write_chr "$here/both" 10 records_10 both

echo "rebuilt $here/{csi,tbi,noindex,endcol,both}"
