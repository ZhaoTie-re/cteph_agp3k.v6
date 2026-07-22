#!/bin/bash
# ---------------------------------------------------------------------------
# Aligned depth and read length of a CRAM, sampled over probe regions.
#
#   $1  CHUNK    TSV of "sample_id <TAB> cram_path", one sample per line
#   $2  REGIONS  chr:start-end[,chr:start-end...]
#   $3  FASTA    reference the CRAMs were encoded against
#   $4  THREADS  parallel samples
#   $5  OUT      output TSV
#
# Output: SAMPLE REGION LEN DEPTH_NODUP DEPTH_WITHDUP READLEN_MEAN READLEN_MIN READLEN_MAX
#
# Why regions rather than the whole CRAM: the audit compares the sheet's
# genome-wide Observed_Depth against the reads, and that comparison only needs a
# read DENSITY, which a sample of the genome estimates fine. `samtools stats` over
# a whole CRAM costs ~1.5 min per sample; a few Mb of probe regions costs seconds.
#
#   DEPTH_WITHDUP  duplicates KEPT. This is the one to compare against
#                  Observed_Depth: that number is derived from the FASTQ, and the
#                  FASTQ still contains the duplicate reads — nothing has marked
#                  them yet. Comparing against the duplicate-free depth would
#                  charge each platform for its duplicate rate on top of the
#                  effect under test.
#   DEPTH_NODUP    samtools coverage at its defaults (drops UNMAP/SECONDARY/
#                  QCFAIL/DUP). Reported so the duplicate burden is visible.
#   READLEN_*      mean/min/max SEQ length over the first READ_SAMPLE reads of the
#                  region. min == max proves a fixed-length library, so the mean
#                  is not an average over a mixture.
# ---------------------------------------------------------------------------
set -euo pipefail

CHUNK=$1; REGIONS=$2; FASTA=$3; THREADS=$4; OUT=$5
READ_SAMPLE=${READ_SAMPLE:-1000}

export FASTA REGIONS READ_SAMPLE

one_sample() {
    id=$1; cram=$2
    for reg in ${REGIONS//,/ }; do
        # column 7 = meandepth = summed per-base depth / region length
        nodup=$(samtools coverage -r "$reg" --reference "$FASTA" "$cram" 2>/dev/null \
                | awk 'NR==2{print $7}')
        withdup=$(samtools coverage -r "$reg" --reference "$FASTA" \
                    --ff UNMAP,SECONDARY,QCFAIL "$cram" 2>/dev/null \
                | awk 'NR==2{print $7}')
        # -F 0x900 drops secondary and supplementary records. Their SEQ is hard-
        # clipped, so they are 30-35 bp fragments of a read rather than reads: left
        # in, they drag the mean down and make every fixed-length library look
        # variable-length.
        rl=$(samtools view -F 0x900 --reference "$FASTA" "$cram" "$reg" 2>/dev/null \
                | head -n "$READ_SAMPLE" \
                | awk '{l=length($10); s+=l; n++; if(l>mx)mx=l; if(mn==""||l<mn)mn=l}
                       END{ if(n) printf "%.1f\t%d\t%d", s/n, mn, mx; else printf "NA\tNA\tNA" }')
        len=$(awk -v r="$reg" 'BEGIN{split(r,a,":"); split(a[2],b,"-"); print b[2]-b[1]+1}')
        [ -n "${nodup:-}" ] && printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$id" "$reg" "$len" "$nodup" "${withdup:-NA}" "${rl:-NA	NA	NA}"
    done
}
export -f one_sample

printf 'SAMPLE\tREGION\tLEN\tDEPTH_NODUP\tDEPTH_WITHDUP\tREADLEN_MEAN\tREADLEN_MIN\tREADLEN_MAX\n' > "$OUT"
# shellcheck disable=SC2016
xargs -a "$CHUNK" -P "$THREADS" -L1 bash -c 'one_sample "$0" "$1"' >> "$OUT"

echo "[cram_region_stats] $(( $(wc -l < "$OUT") - 1 )) rows from $(wc -l < "$CHUNK") samples" >&2
