#!/bin/bash
# ---------------------------------------------------------------------------
# SNP-based association for ONE arm, mirroring analysis/assoc_saige.
#
#   $1  ARM        baseline | rep1 | rep2 | rep3
#   $2  BED_PREFIX PLINK prefix to test
#   $3  NULL       SAIGE null-model prefix (.rda + .varianceRatio.txt)
#   $4  CHROM      chromosome number for step 2
#   $5  OUT        output file
#
# The null model is REUSED, not refitted. It is fit on genome-wide LD-pruned
# markers (5.15M variants); the swap touches 970 of them, so refitting would move
# the null for reasons that have nothing to do with the variants under test, and
# the arms would no longer be comparable. Holding it fixed leaves the genotypes as
# the only thing that differs between arms — which is the entire question.
#
# --AlleleOrder=alt-first. The reference pipeline reads BGEN and passes
# 'ref-first'; this reads PLINK, where the SAME setting silently reverses the sign
# of every BETA. Verified against the reference result at chr16:53887925:T:C:
#     BGEN  ref-first  -> BETA +0.551078
#     PLINK alt-first  -> BETA +0.551078   (matches)
#     PLINK ref-first  -> BETA -0.551078   (reversed)
# Copying the reference's flag across the format change would have inverted every
# effect direction while producing an otherwise entirely normal-looking result.
#
# Assumes PATH/conda are set by the caller (SAIGE).
# ---------------------------------------------------------------------------
set -euo pipefail

ARM=$1; BED_PREFIX=$2; NULL=$3; CHROM=$4; OUT=$5

step2_SPAtests.R \
    --bedFile="${BED_PREFIX}.bed" \
    --bimFile="${BED_PREFIX}.bim" \
    --famFile="${BED_PREFIX}.fam" \
    --AlleleOrder=alt-first \
    --SAIGEOutputFile="${OUT}" \
    --chrom="${CHROM}" \
    --GMMATmodelFile="${NULL}.rda" \
    --varianceRatioFile="${NULL}.varianceRatio.txt" \
    --is_Firth_beta=FALSE \
    --LOCO=TRUE \
    --is_output_moreDetails=TRUE

n=$(( $(wc -l < "${OUT}") - 1 ))
echo "[snp_assoc] ${ARM}: ${n} variants tested" >&2
if [ "${n}" -le 0 ]; then
    echo "[snp_assoc] ${ARM}: SAIGE wrote no result" >&2
    exit 1
fi
