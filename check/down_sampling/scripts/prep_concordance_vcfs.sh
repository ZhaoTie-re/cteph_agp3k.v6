#!/bin/bash
# ---------------------------------------------------------------------------
# Cut ONE replicate down to the refined_core comparison matrix, before and after
# genotype QC.
#
#   $1  REP        replicate label, e.g. rep1
#   $2  NORM_VCF   normalized replicate VCF (04_norm_vcf)
#   $3  VAR_IDS    matrix variant IDs, one per line (00_matrix)
#   $4  SMP_IDS    matrix sample IIDs, one per line (00_matrix)
#   $5  BACKFILL   path to backfill_rgq_to_gq.py
#   $6  FASTA      reference fasta
#   $7  GATK_SIF   GATK singularity image
#   $8  JAVA_OPTS  gatk --java-options string
#
# Outputs:
#   ${REP}.pre_qc.vcf.gz   the matrix, BEFORE genotype QC — the raw down-sampled
#                          calls as the joint caller emitted them
#   ${REP}.post_qc.vcf.gz  the same after adding FORMAT/AF (AlleleFraction),
#                          unphasing, then FILTER_GENOTYPE-style genotype QC with
#                          filtered genotypes set to no-call
#
# Both axes are cut here, so the published VCF *is* the comparison matrix rather
# than something the next step has to re-derive. The matrix itself is defined once
# by build_matrix.sh — this script never touches the .bim.
#
# Mirrors select.auto.par.v6.nf ANNOTATE_AF_NORM_GT + FILTER_GENOTYPE.
# Assumes PATH/conda are already set by the caller (bcftools, bgzip, singularity).
# ---------------------------------------------------------------------------
set -euo pipefail

REP=$1; NORM_VCF=$2; VAR_IDS=$3; SMP_IDS=$4; BACKFILL=$5; FASTA=$6; GATK_SIF=$7; JAVA_OPTS=$8

# Only the matrix samples this replicate actually carries: a replicate holds just
# the down-sampled platforms, so it is a subset of the matrix by design. bcftools
# -S fails on a sample it cannot find, so intersect first rather than assume.
bcftools query -l "${NORM_VCF}" | sort -u > rep_samples.txt
sort -u "${SMP_IDS}" > matrix_samples.sorted.txt
comm -12 rep_samples.txt matrix_samples.sorted.txt > keep_samples.txt
n_keep=$(wc -l < keep_samples.txt)
if [ "${n_keep}" -eq 0 ]; then
    echo "[prep_concordance_vcfs] ${REP}: no sample is in both the replicate and the matrix" >&2
    exit 1
fi
echo "[prep_concordance_vcfs] ${REP}: $(wc -l < rep_samples.txt) called, ${n_keep} on the matrix" >&2

# 1) pre-QC = the matrix, straight from the joint call
bcftools view -i "ID=@${VAR_IDS}" -S keep_samples.txt "${NORM_VCF}" \
    -Oz -o "${REP}.pre_qc.vcf.gz"
bcftools index -f -t "${REP}.pre_qc.vcf.gz"

# 2) Fill GQ from RGQ before any of the QC steps look at it.
#    Forced sites where no sample carries the ALT come out as reference blocks:
#    FORMAT is 'GT:AD:DP:RGQ' with no GQ field at all. The cohort's genotype QC
#    filters on `GQ < 20`, and GATK's JEXL reads a missing GQ as 0, so every one of
#    those confident hom-refs would be no-called for lacking a field rather than
#    for being poor — ~41% of all genotypes here, and only on the down-sampled
#    (i.e. case) samples. See backfill_rgq_to_gq.py. No genotype is changed.
#    pre_qc above stays untouched: it is the raw joint call, and it runs no QC.
bcftools view "${REP}.pre_qc.vcf.gz" | python3 "${BACKFILL}" | bgzip > "${REP}.gq.vcf.gz"
bcftools index -f -t "${REP}.gq.vcf.gz"

# 3) add FORMAT/AF (AlleleFraction) and unphase all genotypes (ANNOTATE_AF_NORM_GT)
singularity exec --bind /LARGE0:/LARGE0 --bind /LARGE1:/LARGE1 \
    "${GATK_SIF}" gatk --java-options "${JAVA_OPTS}" VariantAnnotator \
    -R "${FASTA}" -V "${REP}.gq.vcf.gz" -O "${REP}.af.tmp.vcf.gz" \
    -A AlleleFraction --create-output-variant-index true
bcftools +setGT "${REP}.af.tmp.vcf.gz" -Ou -- -t a -n u \
    | bcftools view -Oz -o "${REP}.af.vcf.gz"
bcftools index -f -t "${REP}.af.vcf.gz"

# 4) genotype QC (FILTER_GENOTYPE): nan->NaN, then set filtered genotypes to no-call
bcftools view "${REP}.af.vcf.gz" | sed 's/nan/NaN/g' | bgzip > "${REP}.af.nan.vcf.gz"
bcftools index -f -t "${REP}.af.nan.vcf.gz"
singularity exec --bind /LARGE0:/LARGE0 --bind /LARGE1:/LARGE1 \
    "${GATK_SIF}" gatk --java-options "${JAVA_OPTS}" VariantFiltration \
    -R "${FASTA}" -V "${REP}.af.nan.vcf.gz" -O "${REP}.post_qc.vcf.gz" \
    --genotype-filter-name "LowGQ"       --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "LowDP"       --genotype-filter-expression "DP < 8" \
    --genotype-filter-name "ABB_outlier" --genotype-filter-expression "isHet == 1 && (AF < 0.2 || AF > 0.8)" \
    --genotype-filter-name "ABB_NaN"     --genotype-filter-expression "AF == 'NaN'" \
    --set-filtered-genotype-to-no-call true \
    --create-output-variant-index true
bcftools index -f -t "${REP}.post_qc.vcf.gz"

rm -f "${REP}.gq.vcf.gz"* "${REP}.af.tmp.vcf.gz"* "${REP}.af.vcf.gz"* "${REP}.af.nan.vcf.gz"*
