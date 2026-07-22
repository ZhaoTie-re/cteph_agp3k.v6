#!/bin/bash
# ---------------------------------------------------------------------------
# Aligned read depth (and read length) of a set of regions, measured straight
# from the CRAMs.
#
#   $1  CHUNK    TSV of "sample_id <TAB> cram_path", one sample per line
#   $2  REGIONS  chr:start-end[,chr:start-end...]
#   $3  FASTA    reference the CRAMs were encoded against
#   $4  THREADS  parallel samples
#   $5  OUT      output TSV
#
# Output TSV: SAMPLE  REGION  LEN  DEPTH_NODUP  DEPTH_WITHDUP  READLEN
#   DEPTH_NODUP    samtools coverage at its defaults, which drop
#                  UNMAP/SECONDARY/QCFAIL/DUP and apply NO mapping- or
#                  base-quality threshold  -> aligned, non-duplicate depth
#   DEPTH_WITHDUP  same but duplicates KEPT -> the pair gives the duplicate burden
#   READLEN        mean SEQ length over the first READ_SAMPLE reads of the region
#
# DEPTH_NODUP is the pool `samtools view -s` subsamples, which is why the
# down-sampling fraction is derived from it rather than from the sample sheet's
# genome-wide Observed_Depth. READLEN is reported per region only to show that it
# is a property of the library, not of where you look.
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
        # mean SEQ length of the first N reads; these libraries are fixed-length
        # (samtools stats reports average == maximum), so a sample is enough
        readlen=$(samtools view --reference "$FASTA" "$cram" "$reg" 2>/dev/null \
                | head -n "$READ_SAMPLE" \
                | awk '{s+=length($10); n++} END{ if(n) printf "%.1f", s/n; else print "NA" }')
        len=$(awk -v r="$reg" 'BEGIN{split(r,a,":"); split(a[2],b,"-"); print b[2]-b[1]+1}')
        [ -n "${nodup:-}" ] && printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$id" "$reg" "$len" "$nodup" "${withdup:-NA}" "${readlen:-NA}"
    done
}
export -f one_sample

printf 'SAMPLE\tREGION\tLEN\tDEPTH_NODUP\tDEPTH_WITHDUP\tREADLEN\n' > "$OUT"
# shellcheck disable=SC2016
xargs -a "$CHUNK" -P "$THREADS" -L1 bash -c 'one_sample "$0" "$1"' >> "$OUT"

echo "[cram_region_depth] $(( $(wc -l < "$OUT") - 1 )) rows from $(wc -l < "$CHUNK") samples" >&2
