#!/usr/bin/env bash
# Step 5: rare driver variants typed directly on one array. The merged array set keeps only sites on both arrays,
# so these come from the per-array sets after sample QC (before merging). WGS genotypes from the full-mainland
# unrelated set; counts are cross-tabulated by 05_rare_direct.py (no sample IDs written).
set -euo pipefail
HERE=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck
MI=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/full_mainland
WGS=$MI/genotype/fixed_model/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready
SQC=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/array/results/02.sample_qc
PL2=/home/b/b37974/plink2
T=${SLURM_CPUS_PER_TASK:-4}
O=$HERE/results/05.rare_direct; mkdir -p $O
printf '%s\n' chr12:112890858:C:T chr11:59576952:G:A chr9:5073770:G:T chr9:5073770:G:A > $O/ids.txt
# WGS: the analysis IDs at the same positions (any alleles), found by position
awk 'BEGIN{OFS="\t"} ($1==12 && $4==112890858) || ($1==11 && $4==59576952) || ($1==9 && $4==5073770) {print $2}' \
    $WGS.bim > $O/wgs.ids
$PL2 --bfile $WGS --extract $O/wgs.ids --export A --out $O/wgs --threads $T
for a in illuminaASA24v1 illuminaJSA24v1; do
  $PL2 --bfile $SQC/$a.acgt.ref.id.selected.sqc --extract $O/ids.txt --export A --out $O/$a --threads $T || true
done
/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 -I $HERE/scripts/05_rare_direct.py
