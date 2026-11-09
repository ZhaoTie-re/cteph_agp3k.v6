#!/usr/bin/env bash
# One joint-call chunk -> the fixed sites with the cohort's own genotype QC.
#
# The steps and expressions are the cohort's (wgs.auto.par/select.auto.par.v6.nf):
#   PREPARE_VCF          norm -m -any --check-ref s, drop ALT="*", ID = CHROM:POS:REF:ALT
#   ANNOTATE_AF_NORM_GT  VariantAnnotator -A AlleleFraction (gatk_latest.sif), +setGT -t a -n u
#   FILTER_GENOTYPE      nan -> NaN, VariantFiltration GQ < 20 | DP < 8 |
#                        het AF outside 0.2-0.8 | AF NaN, filtered genotypes set to no-call
# Between the first two, expand_fixed.py writes exactly one biallelic record per fixed site.
# Site filters (PASS, VQSLOD, MQ) are not re-applied: the fixed sites passed them in the cohort.
#
# usage: genotype_qc.sh <variants.vcf.gz> <allsites.vcf.gz> <fixed_sites.tsv> <chrom> <out prefix>
#        env: PYTHON BCFTOOLS BGZIP TABIX FASTA GATK_SIF EXPAND
set -euo pipefail
var=$1; alls=$2; fixed=$3; chrom=$4; out=$5

$BCFTOOLS norm --multiallelics -any --fasta-ref "$FASTA" --check-ref s "$var" -Ou \
  | $BCFTOOLS filter -e 'ALT="*"' -Ou \
  | $BCFTOOLS annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$out.norm.vcf.gz"
$BCFTOOLS index -f -t "$out.norm.vcf.gz"
$BCFTOOLS index -f -t "$alls"

$PYTHON "$EXPAND" --variants "$out.norm.vcf.gz" --allsites "$alls" --fixed "$fixed" --chrom "$chrom" \
    --out "$out.fixed.vcf"
$BGZIP -f "$out.fixed.vcf"
$TABIX -f -p vcf "$out.fixed.vcf.gz"

singularity exec --bind /LARGE0:/LARGE0 --bind /LARGE1:/LARGE1 "$GATK_SIF" \
    gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4" VariantAnnotator \
    -R "$FASTA" -V "$out.fixed.vcf.gz" -O "$out.af.tmp.vcf.gz" -A AlleleFraction --create-output-variant-index true
$BCFTOOLS +setGT "$out.af.tmp.vcf.gz" -Ou -- -t a -n u | $BCFTOOLS view -Oz -o "$out.af.vcf.gz"
$BCFTOOLS index -f -t "$out.af.vcf.gz"
rm -f "$out.af.tmp.vcf.gz" "$out.af.tmp.vcf.gz.tbi"

$BCFTOOLS view "$out.af.vcf.gz" | sed 's/nan/NaN/g' | $BGZIP > "$out.gtnorm.tmp.vcf.gz"
$BCFTOOLS index -f -t "$out.gtnorm.tmp.vcf.gz"
singularity exec --bind /LARGE0:/LARGE0 --bind /LARGE1:/LARGE1 "$GATK_SIF" \
    gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4" VariantFiltration \
    -R "$FASTA" -V "$out.gtnorm.tmp.vcf.gz" -O "$out.qc.vcf.gz" \
    --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "LowDP" --genotype-filter-expression "DP < 8" \
    --genotype-filter-name "ABB_outlier" --genotype-filter-expression "isHet == 1 && (AF < 0.2 || AF > 0.8)" \
    --genotype-filter-name "ABB_NaN" --genotype-filter-expression "AF == 'NaN'" \
    --set-filtered-genotype-to-no-call true --create-output-variant-index true
rm -f "$out.gtnorm.tmp.vcf.gz" "$out.gtnorm.tmp.vcf.gz.tbi"

n_fixed=$(awk -F'\t' -v c="$chrom" 'NR > 1 && $2 == c' "$fixed" | wc -l)
n_out=$($BCFTOOLS view -H "$out.qc.vcf.gz" | wc -l)
if [ "$n_fixed" -ne "$n_out" ]; then
    echo "ABORT: $chrom has $n_fixed fixed sites but $n_out QC'd records" >&2
    exit 1
fi
echo "[genotype_qc] $chrom: $n_out fixed sites"
