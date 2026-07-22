#!/bin/bash
# ---------------------------------------------------------------------------
# Build ONE arm's gene-based test input: the reference pipeline's own rvtest VCF,
# with this replicate's genotypes swapped in for the down-sampled samples.
#
#   $1  ARM       baseline | rep1 | rep2 | rep3
#   $2  MH_VCF    the reference pipeline's impact-filtered VCF (its rvtest input)
#   $3  REGION    chr:start-end to cut to, in the VCF's own chromosome spelling
#   $4  GENO      'NONE' for the baseline, else the arm's PLINK prefix (08_assoc_genotypes)
#   $5  DS_LIST   'NONE' for the baseline, else that arm's down-sampled sample IDs
#   $6  OUT       output VCF (bgzipped, indexed)
#
# Starting from the reference's own VCF, rather than rebuilding one, is what makes
# the arms comparable: the baseline arm IS the reference input, so it reproduces
# the published p-value exactly, and every replicate differs from it by nothing but
# the genotypes of the down-sampled samples.
#
# Cutting to REGION is safe because rvtest scores each gene from its own variants
# alone; a gene whose variants all lie inside the region gets the same test it
# would in the genome-wide run. Verified: the region cut reproduces the published
# STBD1 result to every digit (2162 / 3 / 3 / 14 / 4.10017e-07).
#
# The variant set is deliberately FIXED at the reference's. Impact annotation does
# not depend on genotypes, so the strata cannot move; MAC could, but letting the
# variant set drift would confound "the genotypes changed" with "a different test
# was run". Per-variant AC is reported so any such drift stays visible.
#
# Assumes PATH/conda are set by the caller (bcftools, plink2).
# ---------------------------------------------------------------------------
set -euo pipefail

ARM=$1; MH_VCF=$2; REGION=$3; GENO=$4; DS_LIST=$5; OUT=$6

# ── the reference's own input, cut to the region ──────────────────────────
bcftools view -r "${REGION}" "${MH_VCF}" -Oz -o ref_region.vcf.gz
bcftools index -f -t ref_region.vcf.gz
n_var=$(bcftools view -H ref_region.vcf.gz | wc -l)
n_smp=$(bcftools query -l ref_region.vcf.gz | wc -l)
if [ "${n_var}" -eq 0 ]; then
    echo "[gene_assoc_prep] ${ARM}: no variant of ${MH_VCF} lies in ${REGION}" >&2
    exit 1
fi
echo "[gene_assoc_prep] ${ARM}: reference region carries ${n_var} variants x ${n_smp} samples" >&2

if [ "${GENO}" = "NONE" ]; then
    # The baseline arm is the reference input itself, untouched.
    cp ref_region.vcf.gz "${OUT}"
    bcftools index -f -t "${OUT}"
    echo "[gene_assoc_prep] ${ARM}: baseline — reference genotypes, nothing swapped" >&2
else
    # Swap only samples the reference actually kept: its own QC dropped a few of
    # the down-sampled ones, and bcftools would fail on a name it cannot find.
    bcftools query -l ref_region.vcf.gz | sort -u > ref_samples.txt
    sort -u "${DS_LIST}" > ds_all.txt
    comm -12 ref_samples.txt ds_all.txt > swap_samples.txt
    n_swap=$(wc -l < swap_samples.txt)
    if [ "${n_swap}" -eq 0 ]; then
        echo "[gene_assoc_prep] ${ARM}: none of the down-sampled samples survive in the reference VCF" >&2
        exit 1
    fi
    echo "[gene_assoc_prep] ${ARM}: swapping ${n_swap} of $(wc -l < ds_all.txt) down-sampled samples" >&2

    # The variants to swap at, and the REF each must keep. The ID is
    # CHROM:POS:REF:ALT, so REF is field 3 — forcing it stops the PLINK->VCF export
    # from choosing its own and breaking the merge.
    bcftools view -H ref_region.vcf.gz | cut -f3 | sort -u > vars.txt
    awk -F: '{print $0"\t"$3}' vars.txt > ref_allele.txt
    awk 'NR==FNR{s[$1];next} ($2 in s){print $1"\t"$2}' swap_samples.txt "${GENO}.fam" > swap_keep.txt

    # --output-chr must match the reference VCF's spelling, or bcftools merge sees
    # two different chromosomes, merges nothing, and emits every variant twice.
    chr_style=$(bcftools view -H ref_region.vcf.gz | head -1 | cut -f1 | grep -q '^chr' && echo chrM || echo 26)
    plink2 --bfile "${GENO}" \
        --keep swap_keep.txt --extract vars.txt \
        --ref-allele force ref_allele.txt 2 1 \
        --output-chr "${chr_style}" \
        --export vcf bgz id-paste=iid --out rep_part
    bcftools index -f -t rep_part.vcf.gz

    # everyone the reference keeps, minus the swapped ones
    bcftools view -S ^swap_samples.txt --force-samples ref_region.vcf.gz -Oz -o ref_rest.vcf.gz
    bcftools index -f -t ref_rest.vcf.gz

    bcftools merge ref_rest.vcf.gz rep_part.vcf.gz -Oz -o "${OUT}"
    bcftools index -f -t "${OUT}"

    n_out_var=$(bcftools view -H "${OUT}" | wc -l)
    n_out_smp=$(bcftools query -l "${OUT}" | wc -l)
    if [ "${n_out_var}" -ne "${n_var}" ] || [ "${n_out_smp}" -ne "${n_smp}" ]; then
        echo "[gene_assoc_prep] ${ARM}: got ${n_out_var} variants x ${n_out_smp} samples," \
             "expected ${n_var} x ${n_smp}. A variant count that is a multiple of the" \
             "expectation means the merge matched nothing and kept both copies." >&2
        exit 1
    fi
fi

# Per-variant allele counts, so a genotype change that moves MAC is visible rather
# than buried inside a gene-level p-value.
bcftools +fill-tags "${OUT}" -Ou -- -t AC,AN 2>/dev/null \
    | bcftools query -f "${ARM}\t%ID\t%AC\t%AN\n" > "${ARM}.ac.tsv"
echo "[gene_assoc_prep] ${ARM}: wrote ${OUT}" >&2
