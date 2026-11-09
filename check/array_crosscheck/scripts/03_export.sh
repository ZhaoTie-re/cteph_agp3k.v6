#!/usr/bin/env bash
# Step 3: genotypes around the signals (WGS: leads, credible-set members, shared array sites; array: shared sites
# and the MHC), the WGS association at the leads in the same people, and concordance by case status.
set -euo pipefail

HERE=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck
MI=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/full_mainland
WGS=$MI/genotype/fixed_model/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready
# the array set after QC (sample call rate, merge, biallelic, variant missingness, autosomes, MAF, HWE) and
# before the ToMMo filter: 1,750 samples, 409,088 variants
ARR=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/array/results/10.variant_qc_filtered/cteph_agp3k.array.sqc.vqc
PHENO=$MI/phenotype/pheno.tsv
COV=$MI/covariates/bbj_mainland_pc.sex.tsv
PL2=/home/b/b37974/plink2
PL1=/home/b/b37974/plink
T=${SLURM_CPUS_PER_TASK:-4}
O=$HERE/results/01.inputs; A=$HERE/results/02.assoc; C=$HERE/results/03.concordance

$PL2 --bfile $WGS --keep $O/keep.txt --extract $O/wgs_window.ids --export A --out $O/wgs.window --threads $T
$PL2 --bfile $ARR --keep $O/keep.txt --extract $O/arr_window.ids --export A --out $O/arr.window --threads $T
$PL2 --bfile $ARR --keep $O/keep.txt --extract range $O/mhc.range --maf 0.01 --export A --out $O/arr.mhc --threads $T
$PL2 --bfile $WGS --keep $O/keep.txt --extract $O/wgs_window.ids --pheno $PHENO --pheno-name PHENO1 --1 \
     --covar $COV --covar-name SEX PC1_AVG-PC10_AVG --covar-variance-standardize \
     --glm firth-fallback hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p \
     --out $A/wgs.window --threads $T
grep -q "11 covariates loaded" $A/wgs.window.log
# the array on its own near every lead (also sites the WGS set lacks, e.g. leads that failed the WGS variant QC)
$PL2 --bfile $ARR --keep $O/keep.txt --extract range $O/lead_windows.range --pheno $PHENO --pheno-name PHENO1 --1 \
     --covar $COV --covar-name SEX PC1_AVG-PC10_AVG --covar-variance-standardize \
     --glm firth-fallback hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p \
     --out $A/arr.window --threads $T
grep -q "11 covariates loaded" $A/arr.window.log

# concordance of non-missing calls, cases and controls apart (the log gives the totals)
# (plink 1.9 cannot filter samples while it merges, so each group is subset first)
for G in cases controls; do
  for S in wgs arr; do
    $PL2 --bfile $O/$S.shared --keep $O/keep.$G.txt --make-bed --out $C/tmp.$S.$G --threads $T
  done
  $PL1 --bfile $C/tmp.wgs.$G --bmerge $C/tmp.arr.$G --merge-mode 7 --out $C/wgs_vs_arr.$G --threads $T || true
  grep -i "concordan" $C/wgs_vs_arr.$G.log
  rm -f $C/tmp.*.$G.*
done
