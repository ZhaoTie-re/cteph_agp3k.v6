#!/bin/bash
# ---------------------------------------------------------------------------
# Gene-based test for ONE arm x ONE method, mirroring analysis/assoc_rvtest.
#
#   $1  ARM        baseline | rep1 | rep2 | rep3
#   $2  METHOD     tag, e.g. skato / cmc / zeggini
#   $3  METHOD_OPT rvtest option, e.g. '--kernel skato'
#   $4  VCF        this arm's test VCF (gene_assoc_prep.sh)
#   $5  PHENO      rvtest-format phenotype file
#   $6  PHENO_NAME phenotype column
#   $7  COVAR      rvtest-format covariate file
#   $8  COVAR_NAME covariate columns, comma-separated
#   $9  GENEFILE   refFlat, in the VCF's chromosome spelling
#   $10 OUT_PREFIX output prefix
#
# The phenotype/covariate/refFlat files are the reference pipeline's OWN prepared
# copies, not the upstream originals: rvtest wants 'fid iid fatid matid sex <pheno>'
# and a refFlat whose chromosomes match the VCF's, and the raw files are neither.
# Handing it the originals fails loudly on the covariates but merely *warns* on the
# phenotype — it drops every sample and writes an empty result that looks like a
# test that found nothing.
#
# Assumes PATH is set by the caller (rvtest).
# ---------------------------------------------------------------------------
set -euo pipefail

ARM=$1; METHOD=$2; METHOD_OPT=$3; VCF=$4
PHENO=$5; PHENO_NAME=$6; COVAR=$7; COVAR_NAME=$8; GENEFILE=$9; OUT_PREFIX=${10}

rvtest \
    --inVcf "${VCF}" \
    --pheno "${PHENO}" --pheno-name "${PHENO_NAME}" \
    --covar "${COVAR}" --covar-name "${COVAR_NAME}" \
    --geneFile "${GENEFILE}" \
    --out "${OUT_PREFIX}" \
    --noweb \
    --numThread 4 \
    ${METHOD_OPT}

# rvtest reports "dropped due to missing phenotype" as a warning and still exits 0,
# leaving an empty .assoc behind. Refuse that instead of reporting it as a result.
n=$(cat "${OUT_PREFIX}".*.assoc 2>/dev/null | tail -n +2 | wc -l)
echo "[gene_assoc_run] ${ARM}/${METHOD}: ${n} gene(s) tested" >&2
if [ "${n}" -eq 0 ]; then
    echo "[gene_assoc_run] ${ARM}/${METHOD}: rvtest produced no gene result. Tail of its log:" >&2
    tail -20 "${OUT_PREFIX}.log" >&2 || true
    exit 1
fi
