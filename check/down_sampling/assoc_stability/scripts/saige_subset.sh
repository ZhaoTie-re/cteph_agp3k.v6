#!/bin/bash
# ---------------------------------------------------------------------------
# Run the MAIN ANALYSIS on a sample subset: SAIGE step 1 + step 2, exactly as
# analysis/assoc_saige runs it, on whichever samples this subset names.
#
#   $1  TAG        subset label, e.g. minus_DNBSeq-T7
#   $2  KEEP       FID<TAB>IID of the samples to keep
#   $3  BASE       random_model PLINK prefix (step 1's genotypes)
#   $4  PRUNE_IN   LD-pruned marker list from the reference run
#   $5  TEST_BED   PLINK prefix to test in step 2 (the region)
#   $6  MERGE_PY   the reference pipeline's merge_pheno_cov.py
#   $7  PHENO_FILE raw phenotype file
#   $8  COV_FILE   raw covariate file
#   $9  PHENO_COL  phenotype column
#   $10 COVAR      covariate column list, comma-separated
#   $11 CHROM      chromosome for step 2
#   $12 THREADS    step 1 threads
#   $13 OUT        output .assoc.txt
#
# WHY SAIGE AND NOT A CRUDE FISHER TEST
#   Fisher on allele counts would cost seconds instead of half an hour, and would
#   be a different analysis, not a cheaper one:
#     - random_model is NOT relatedness-pruned. It holds 2,193 samples including
#       related pairs, deliberately, because the GRM is what absorbs them. Fisher
#       treats relatives as independent and returns an anticonservative p.
#     - the model adjusts SEX + 10 PCs, and the PCs differ between cases and
#       controls here (PC1 0.0066 vs 0.0079): they come from different collections.
#     - the published effect is exp(BETA) from this model. A crude OR is another
#       estimand — 1.62 against the model's 1.74 at the lead variant — and putting
#       the two side by side invites exactly the comparison that must not be made.
#
# WHY THE NULL IS REFIT PER SUBSET
#   The null carries the fitted values of the samples it saw. Dropping a platform
#   and keeping the full-cohort null answers a question nobody asked; refitting is
#   what "the study without this platform" means. ~26 min per subset, and the
#   subsets run in parallel.
#
# Assumes PATH/conda are set by the caller (SAIGE, plink2, python3).
# ---------------------------------------------------------------------------
set -euo pipefail

TAG=$1; KEEP=$2; BASE=$3; PRUNE_IN=$4; TEST_BED=$5
MERGE_PY=$6; PHENO_FILE=$7; COV_FILE=$8; PHENO_COL=$9; COVAR=${10}
CHROM=${11}; THREADS=${12}; OUT=${13}

echo "[saige_subset] ${TAG}: $(wc -l < "${KEEP}") samples requested" >&2

# ── a case-control model needs both classes ────────────────────────────────
awk 'NR==FNR{k[$2];next} ($2 in k){c[$6]++} END{for (p in c) print p, c[p]}' \
    "${KEEP}" "${BASE}.fam" | sort > pheno_mix.txt
n_class=$(wc -l < pheno_mix.txt)
n_case=$(awk '$1=="2"{print $2}' pheno_mix.txt)
echo "[saige_subset] ${TAG}: phenotype mix (1=control 2=case): $(tr '\n' ' ' < pheno_mix.txt)" >&2
if [ "${n_class}" -lt 2 ]; then
    echo "[saige_subset] ${TAG}: only one phenotype class present — a case-control model" >&2
    echo "  cannot be fit. In this cohort that is what 'within one platform' means, since" >&2
    echo "  no platform holds both cases and controls." >&2
    exit 1
fi
# Not fatal, but the reader has to know: 11 covariates on a handful of cases gives
# an interval so wide it says nothing, and that is a real result, not a failure.
if [ "${n_case:-0}" -lt 50 ]; then
    echo "[saige_subset] ${TAG}: WARNING — only ${n_case} cases for 11 covariates." >&2
    echo "  Expect a wide confidence interval; read it as imprecision, not as absence." >&2
fi

# ── step 1: the null, refit on THIS subset ─────────────────────────────────
plink2 --bfile "${BASE}" --keep "${KEEP}" --extract "${PRUNE_IN}" \
    --make-bed --out step1_geno --threads "${THREADS}"

# The reference pipeline's own merge script, so the covariate table is built the
# same way rather than a lookalike of it.
python3 "${MERGE_PY}" \
    --pheno "${PHENO_FILE}" \
    --cov "${COV_FILE}" \
    --pheno_col "${PHENO_COL}" \
    --cov_list "${COVAR}" \
    --sex_col "SEX" \
    --out merged_pheno_cov.txt

step1_fitNULLGLMM.R \
    --plinkFile=step1_geno \
    --phenoFile=merged_pheno_cov.txt \
    --phenoCol="${PHENO_COL}" \
    --covarColList="${COVAR}" \
    --sexCol=SEX \
    --sampleIDColinphenoFile=IID \
    --traitType=binary \
    --outputPrefix="${TAG}.null" \
    --nThreads="${THREADS}" \
    --isDiagofKinSetAsOne=True \
    --numRandomMarkerforVarianceRatio=200 \
    --skipVarianceRatioEstimation=FALSE \
    --useSparseGRMtoFitNULL=FALSE \
    --IsOverwriteVarianceRatioFile=TRUE \
    --isCovariateOffset=FALSE

# ── step 2: the region, on the same subset ────────────────────────────────
plink2 --bfile "${TEST_BED}" --keep "${KEEP}" --make-bed --out step2_geno

# --AlleleOrder=alt-first for PLINK input. The reference reads BGEN and passes
# 'ref-first'; the same flag on PLINK reverses the sign of every BETA. Verified
# against the reference at chr16:53887925:T:C — alt-first reproduces +0.551078,
# ref-first gives -0.551078.
step2_SPAtests.R \
    --bedFile=step2_geno.bed \
    --bimFile=step2_geno.bim \
    --famFile=step2_geno.fam \
    --AlleleOrder=alt-first \
    --SAIGEOutputFile="${OUT}" \
    --chrom="${CHROM}" \
    --GMMATmodelFile="${TAG}.null.rda" \
    --varianceRatioFile="${TAG}.null.varianceRatio.txt" \
    --is_Firth_beta=FALSE \
    --LOCO=TRUE \
    --is_output_moreDetails=TRUE

n=$(( $(wc -l < "${OUT}") - 1 ))
echo "[saige_subset] ${TAG}: ${n} variants tested" >&2
[ "${n}" -gt 0 ] || { echo "[saige_subset] ${TAG}: SAIGE wrote no result" >&2; exit 1; }
