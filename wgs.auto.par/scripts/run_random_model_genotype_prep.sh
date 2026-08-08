#!/bin/zsh
# =============================================================================
# run_random_model_genotype_prep.sh
# -----------------------------------------------------------------------------
# Purpose : Random-model genotype: all samples x fixed-model common (MAF>=thr) variants (GRM).
# Project : cteph_agp3k.v6 WGS pipeline (wgs.auto.par/select.auto.par.v6.nf)
# Used by : PREPARE_RANDOM_MODEL_GENOTYPE
# =============================================================================
# Steps:
# 1) Extract variants with MAF >= threshold (from fixed-model variant list)
#    from PopGMM-subset genotype for random-effects model.
# =============================================================================
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <popgmm_bfile_prefix> <maf_ge_variants_list> <out_prefix> <threads>" >&2
  exit 1
fi

popgmm_bfile_prefix="$1"
maf_ge_variants_list="$2"
out_prefix="$3"
threads="$4"

log_file="${out_prefix}.random_model_prep.log.txt"

if [[ ! -f "${popgmm_bfile_prefix}.bed" ]]; then
  echo "[ERROR] Input bfile not found: ${popgmm_bfile_prefix}.bed" >&2
  exit 1
fi

if [[ ! -f "${maf_ge_variants_list}" ]]; then
  echo "[ERROR] MAF >= threshold variants list not found: ${maf_ge_variants_list}" >&2
  exit 1
fi

# Extract variants with MAF >= threshold for random model
plink2 \
  --bfile "${popgmm_bfile_prefix}" \
  --extract "${maf_ge_variants_list}" \
  --make-bed \
  --out "${out_prefix}" \
  --threads "${threads}"

# Generate summary log
plink_version=$(plink2 --version 2>/dev/null | head -1 || echo "plink2 (version unavailable)")
pop_n=$(wc -l < "${popgmm_bfile_prefix}.fam")
pop_v=$(wc -l < "${popgmm_bfile_prefix}.bim")
random_n=$(wc -l < "${out_prefix}.fam")
random_v=$(wc -l < "${out_prefix}.bim")
extracted_v=$(wc -l < "${maf_ge_variants_list}")

{
  echo "=============================================================================="
  echo "Random-model genotype preparation"
  echo "=============================================================================="
  echo "TIMESTAMP        : $(date '+%Y-%m-%d %H:%M:%S')"
  echo "TOOL             : ${plink_version}"
  echo "INPUT_BFILE      : ${popgmm_bfile_prefix}  (${pop_n} samples / ${pop_v} variants)"
  echo "VARIANT_LIST_SRC : ${maf_ge_variants_list}  (fixed-model common part, ${extracted_v} variants)"
  echo ""
  echo "Step 1 | NO sample removal — all ${pop_n} set samples retained"
  echo "Step 2 | Extract common (MAF>=threshold) variants (plink2 --extract)"
  echo "  variants: ${pop_v} -> ${random_v}"
  echo "  output bfile                     : ${out_prefix}.{bed,bim,fam}  (${random_n} / ${random_v})"
  echo ""
  echo "Interpretation:"
  echo "  Random-model genotype = ALL set samples x common (MAF>=threshold) variants."
  echo "  Variant set is identical to the fixed-model common part; used to build the"
  echo "  GRM / random effect (relatedness handled by the model, not by sample removal)."
  echo "=============================================================================="
} > "$log_file"

echo "[OK] Random model genotype prepared: $out_prefix"
