#!/bin/bash
# ---------------------------------------------------------------------------
# Genotype concordance of all replicates (pre_qc & post_qc) against the cohort's
# own genotypes, on the refined_core comparison matrix. Produces ONE combined log
# and two figures via genotype_concordance.py.
#
#   $1  REF      refined_core PLINK prefix (.bed/.bim/.fam) — the truth genotypes
#   $2  VAR_IDS  matrix variant IDs, one per line (00_matrix)
#   $3  SMP_IDS  matrix sample IIDs, one per line (00_matrix)
#   $4  PY       path to genotype_concordance.py
#
# Expects the per-replicate VCFs staged in the CWD:
#   <rep>.pre_qc.vcf.gz  and  <rep>.post_qc.vcf.gz
#
# The matrix is defined once by build_matrix.sh; this script never re-derives it
# from the .bim. It does intersect the samples once more against the replicates,
# because only the down-sampled platforms were called at all.
#
# truth is exported as ALT dosage (--export-allele forces the ALT); test GT come
# straight from the VCFs (biallelic, so GT is ALT-relative) — orientation-safe.
# Assumes PATH/conda are set by the caller (plink2, bcftools, python3).
# ---------------------------------------------------------------------------
set -euo pipefail

REF=$1; VAR_IDS=$2; SMP_IDS=$3; PY=$4

# ── scored samples = matrix ∩ what the replicates actually contain ──────────
# prep_concordance_vcfs.sh already cut the VCFs to the matrix, so this is normally
# a no-op — but deriving it from the VCFs rather than assuming keeps the truth and
# test sides provably identical.
one_vcf=$(ls *.pre_qc.vcf.gz | head -1)
bcftools query -l "${one_vcf}" | sort -u > rep_samples.txt
sort -u "${SMP_IDS}" > matrix_samples.sorted.txt
comm -12 rep_samples.txt matrix_samples.sorted.txt > shared_samples.txt
echo "[genotype_concordance] matrix $(wc -l < matrix_samples.sorted.txt) samples;" \
     "replicates carry $(wc -l < rep_samples.txt); scoring $(wc -l < shared_samples.txt)" >&2

# ── the ALT allele to count, per variant ───────────────────────────────────
# IDs are CHROM:POS:REF:ALT, so ALT is the last colon-separated field. Naming it
# explicitly makes the dosage ALT-relative, immune to which allele PLINK happens
# to treat as A1.
awk -F: '{print $0"\t"$NF}' "${VAR_IDS}" > alt_allele.txt
awk 'NR==FNR{s[$1];next} ($2 in s){print $1,$2}' shared_samples.txt "${REF}.fam" > keep_ref.txt

# TRUTH = the cohort's QC'd genotypes -> ALT dosage 0/1/2/NA
plink2 --bfile "${REF}" \
    --keep keep_ref.txt --extract "${VAR_IDS}" \
    --export A --export-allele alt_allele.txt \
    --out truth

# TEST = each condition's genotypes as a long table (SAMPLE  ID  GT)
for v in *.pre_qc.vcf.gz *.post_qc.vcf.gz; do
    label=$(basename "${v}" .vcf.gz)          # e.g. rep1.pre_qc
    bcftools query -f '[%SAMPLE\t%ID\t%GT\n]' -S shared_samples.txt "${v}" > "${label}.gt.tsv"
done

# confusion matrices + metrics + missingness -> one log + two figures
python3 "${PY}" \
    --truth-raw truth.raw \
    --gt-dir . \
    --out-log genotype_concordance.log \
    --out-fig genotype_concordance_matrices.png \
    --out-fig-missing genotype_missingness.png
