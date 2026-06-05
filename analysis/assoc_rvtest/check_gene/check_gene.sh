#!/bin/bash
# =========================================================================================
#  check_gene.sh  ·  Per-gene deep-dive for the rvtests-rvat pipeline (cteph_agp3k.v6)
# -----------------------------------------------------------------------------------------
#  Locates this gene's rvtest outputs (raw assoc + CMC/SKAT-O/Zeggini FDR + impact VCF)
#  and runs check_gene_detail.py to report per-variant / per-sample details, cumulative
#  MAC/MAF in cases vs controls, and ToMMo allele frequencies.
#
#  Usage: ./check_gene.sh <Gene> [Impact] [Sample_Group]
#     Impact        : high | moderate_high | low_moderate_high     (default: moderate_high)
#     Sample_Group  : both | case | control                        (default: both)
# =========================================================================================

# Always run under bash, even if launched as `zsh check_gene.sh` (other shells may start
# with a stripped-down PATH and different command hashing). Prefer running: ./check_gene.sh
if [ -z "${BASH_VERSION:-}" ]; then exec /bin/bash "$0" "$@"; fi
set -uo pipefail

# Make sure core utilities (and the project tools bcftools/tabix) resolve even under a
# stripped-down non-interactive PATH. plink2 is referenced by absolute path; python comes
# from the conda env activated below.
export PATH="/usr/bin:/bin:/home/b/b37974/bcftools:/home/b/b37974/htslib-1.9:${PATH:-}"
hash -r 2>/dev/null || true

# --- Colors ---
RED='\033[0;31m'; GREEN='\033[0;32m'; BLUE='\033[0;34m'
YELLOW='\033[1;33m'; CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

# --- Arguments ---
GENE_NAME=${1:-}
IMPACT=${2:-moderate_high}
SAMPLE_GROUP=${3:-both}

if [ -z "$GENE_NAME" ]; then
    echo -e "${RED}Usage: $0 <Gene_Name> [Impact] [Sample_Group]${NC}"
    echo "  Impact       : high | moderate_high (default) | low_moderate_high"
    echo "  Sample_Group : both (default) | case | control"
    exit 1
fi

# --- Pipeline layout (new rvtest pipeline) ---
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
ANALYSIS_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/assoc_rvtest"
ROOT_DIR="${ANALYSIS_DIR}/results"
GTQC_DIR="${ROOT_DIR}/00.filter_genotype"
PREP_DIR="${ROOT_DIR}/01.rvtest_prepare"
VCF_DIR="${ROOT_DIR}/03.info_filter"
ASSOC_DIR="${ROOT_DIR}/04.rvtest_run"
POST_DIR="${ROOT_DIR}/05.post_process"

# --- Shared external resources ---
TOMMO_VCF="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
PLINK2_PATH="/home/b/b37974/plink2_alpha6/plink2"
# Fallback genotype matrix if the gtqc-filtered set is unavailable (filterGenotype was off):
INPUT_PLINK="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold"

# --- Normalise the impact tag (new tags carry no 'impact_' prefix) ---
case "$IMPACT" in
    high|impact_high)                              FILTER="high" ;;
    moderate_high|mod_high|impact_moderate_high)   FILTER="moderate_high" ;;
    low_moderate_high|all|impact_low_moderate_high) FILTER="low_moderate_high" ;;
    *) echo -e "${YELLOW}[WARN] Unknown impact '$IMPACT'; using it verbatim.${NC}"; FILTER="$IMPACT" ;;
esac

# --- Locate result files by the stable <filter>/<method>/ layout (prefix-agnostic) ---
ASSOC_FILE=$(ls  "${ASSOC_DIR}/${FILTER}/skato/"*.SkatO.assoc                    2>/dev/null | head -n1)
CMC_FILE=$(ls    "${POST_DIR}/${FILTER}/cmc/"*.CMC.filtered.fdr.assoc           2>/dev/null | head -n1)
SKATO_FILE=$(ls  "${POST_DIR}/${FILTER}/skato/"*.SkatO.filtered.fdr.assoc       2>/dev/null | head -n1)
ZEGGINI_FILE=$(ls "${POST_DIR}/${FILTER}/zeggini/"*.Zeggini.filtered.fdr.assoc  2>/dev/null | head -n1)
VCF_FILE=$(ls    "${VCF_DIR}/"*."${FILTER}".vcf.gz                              2>/dev/null | head -n1)
PHENO_FILE=$(ls  "${PREP_DIR}/"*.pheno_rvt.tsv                                  2>/dev/null | head -n1)
COVAR_FILE=$(ls  "${PREP_DIR}/"*.covar_rvt.tsv                                  2>/dev/null | head -n1)
REFFLAT_FILE="${PREP_DIR}/refFlat.hg38.nochr.txt.gz"
# Covariates used by the pipeline (must match rv_test.nf params.covarName)
COVAR_NAME="sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg"

# Genotype matrix that went into the analysis: prefer the gtqc-filtered PLINK set.
GTQC_BED=$(ls "${GTQC_DIR}/"*.gtqc.bed 2>/dev/null | head -n1)
if [ -n "$GTQC_BED" ]; then
    PLINK_PREFIX="${GTQC_BED%.bed}"
    GENO_SRC="gtqc-filtered (00.filter_genotype)"
else
    PLINK_PREFIX="$INPUT_PLINK"
    GENO_SRC="original input (unfiltered)"
fi

# --- Banner ---
echo -e "${CYAN}==========================================${NC}"
echo -e "${BOLD}Target Gene     : ${GREEN}${GENE_NAME}${NC}"
echo -e "${BOLD}Impact Stratum  : ${BLUE}${FILTER}${NC}"
echo -e "${BOLD}Sample Scope    : ${BLUE}${SAMPLE_GROUP}${NC}"
echo -e "${BOLD}Genotype Source : ${BLUE}${GENO_SRC}${NC}"
echo -e "${CYAN}==========================================${NC}"

# --- Validate critical inputs ---
MISSING=0
for pair in "Raw assoc:$ASSOC_FILE" "Impact VCF:$VCF_FILE" "Pheno:$PHENO_FILE" "PLINK bed:${PLINK_PREFIX}.bed"; do
    label="${pair%%:*}"; path="${pair#*:}"
    if [ -z "$path" ] || [ ! -e "$path" ]; then
        echo -e "${RED}[ERROR] ${label} not found for stratum '${FILTER}'${NC}  ($path)"
        MISSING=1
    fi
done
for pair in "CMC FDR:$CMC_FILE" "SKAT-O FDR:$SKATO_FILE" "Zeggini FDR:$ZEGGINI_FILE" "refFlat:$REFFLAT_FILE"; do
    label="${pair%%:*}"; path="${pair#*:}"
    if [ -z "$path" ] || [ ! -e "$path" ]; then
        echo -e "${YELLOW}[WARN] ${label} not found (gene-level stat will be NA): $path${NC}"
    fi
done
[ "$MISSING" -eq 1 ] && { echo -e "${RED}Aborting: a critical input is missing.${NC}"; exit 1; }

# --- Output locations ---
OUT_BASE="${SCRIPT_DIR}/output/${GENE_NAME}"
mkdir -p "${OUT_BASE}"
TMP_DIR="${OUT_BASE}/tmp/${FILTER}_$$"
mkdir -p "${TMP_DIR}"
LOG_FILE="${OUT_BASE}/${GENE_NAME}.${FILTER}.summary.txt"
rm -f "${LOG_FILE}"

# Activate the conda env robustly: source conda.sh to define the `conda` shell function
# (needed for `conda activate` in a non-interactive shell), then activate.
for csh in /home/b/b37974/anaconda3/etc/profile.d/conda.sh "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    [ -f "$csh" ] && source "$csh" && break
done
conda activate cteph_geno_pro 2>/dev/null || source activate cteph_geno_pro 2>/dev/null

# --- Run the detail analyser ---
python "${SCRIPT_DIR}/check_gene_detail.py" \
    --gene "${GENE_NAME}" \
    --assoc-file "${ASSOC_FILE}" \
    --burden-file "${CMC_FILE}" \
    --skato-file "${SKATO_FILE}" \
    --zeggini-file "${ZEGGINI_FILE}" \
    --vcf-file "${VCF_FILE}" \
    --plink-prefix "${PLINK_PREFIX}" \
    --tommo-vcf "${TOMMO_VCF}" \
    --pheno-file "${PHENO_FILE}" \
    --refflat-file "${REFFLAT_FILE}" \
    --covar-file "${COVAR_FILE}" \
    --covar-name "${COVAR_NAME}" \
    --plink2-path "${PLINK2_PATH}" \
    --out-dir "${TMP_DIR}" \
    --out-log "${LOG_FILE}" \
    --sample-group "${SAMPLE_GROUP}" \
    --group-name "${FILTER}"

rm -rf "${TMP_DIR}"

# --- Final report ---
echo -e ""
echo -e "${CYAN}============================================================${NC}"
echo -e "${BOLD}                   ANALYSIS COMPLETE                        ${NC}"
echo -e "${CYAN}============================================================${NC}"
echo -e "  > Gene Target   : ${GREEN}${GENE_NAME}${NC}"
echo -e "  > Impact Stratum: ${BLUE}${FILTER}${NC}"
echo -e "  > Sample Scope  : ${YELLOW}${SAMPLE_GROUP}${NC}"
echo -e "  > Output Dir    : ${OUT_BASE}"
echo -e "  > Summary Log   : ${GREEN}${LOG_FILE}${NC}"

CASE_DETAILS="${OUT_BASE}/${GENE_NAME}.${FILTER}.case.sample_details.tsv"
CTRL_DETAILS="${OUT_BASE}/${GENE_NAME}.${FILTER}.control.sample_details.tsv"
[ -f "$CASE_DETAILS" ] && echo -e "  > Case Details  : ${GREEN}${CASE_DETAILS}${NC}"
[ -f "$CTRL_DETAILS" ] && echo -e "  > Ctrl Details  : ${GREEN}${CTRL_DETAILS}${NC}"
echo -e "${CYAN}============================================================${NC}"
