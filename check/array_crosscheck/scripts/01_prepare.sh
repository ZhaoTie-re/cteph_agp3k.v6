#!/usr/bin/env bash
# Step 1: the same people on two technologies. Unrelated full-mainland genomes that also have SNP-array genotypes;
# the WGS and array sets restricted to them and to the variant IDs both carry; association, missingness, HWE, and
# genotype concordance on each. Inputs are read only.
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
mkdir -p $O $A $C

# samples: unrelated full mainland (fixed model) with array genotypes; cases and controls separately
awk 'NR==FNR{a[$2];next} ($2 in a){print $1"\t"$2}' $ARR.fam $WGS.fam > $O/keep.txt
awk 'NR==FNR{k[$2];next} FNR>1 && ($2 in k) && $3==1{print $1"\t"$2}' $O/keep.txt $PHENO > $O/keep.cases.txt
awk 'NR==FNR{k[$2];next} FNR>1 && ($2 in k) && $3==0{print $1"\t"$2}' $O/keep.txt $PHENO > $O/keep.controls.txt

# variants: IDs (chr:pos:ref:alt) carried by both
awk 'NR==FNR{a[$2];next} ($2 in a){print $2}' $ARR.bim $WGS.bim > $O/shared.ids

for S in wgs arr; do
  IN=$WGS; [ $S = arr ] && IN=$ARR
  $PL2 --bfile $IN --keep $O/keep.txt --extract $O/shared.ids --make-bed --out $O/$S.shared --threads $T
  # the same model on both: logistic, sex + 10 BBJ-mainland PCs, Firth when the fit separates
  $PL2 --bfile $O/$S.shared --pheno $PHENO --pheno-name PHENO1 --1 \
       --covar $COV --covar-name SEX PC1_AVG-PC10_AVG --covar-variance-standardize \
       --glm firth-fallback hide-covar \
             cols=chrom,pos,ref,alt,a1freq,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p \
       --out $A/$S --threads $T
  for G in cases controls; do
    $PL2 --bfile $O/$S.shared --keep $O/keep.$G.txt --missing variant-only --out $C/$S.$G --threads $T
  done
  $PL2 --bfile $O/$S.shared --keep $O/keep.controls.txt --hardy --out $C/$S.controls --threads $T
done
# the covariates must have entered the model (plink2 does not stop when a named column is absent)
grep -q "11 covariates loaded" $A/wgs.log && grep -q "11 covariates loaded" $A/arr.log

# genotype concordance: every non-missing call that differs between the two technologies
$PL1 --bfile $O/wgs.shared --bmerge $O/arr.shared --merge-mode 7 --out $C/wgs_vs_arr --threads $T || true
test -s $C/wgs_vs_arr.diff
