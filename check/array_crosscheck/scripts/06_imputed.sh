#!/usr/bin/env bash
# Step 6: the imputed arrays (Minimac4, 2025-10-28). The reference panel is AGP3K whole genomes (the user,
# 2026-10-30), i.e. the controls' own sequences: imputed controls largely reproduce their WGS calls. The cases'
# imputed genotypes come from their array and the panel haplotypes only, without their DNBSeq-T7 reads, so they are
# the independent part. Exports, for the same 1,683 people:
#   (a) rare variants (WGS minor allele count 1-30) within +/-500 kb of each rare driver, to measure carrier
#       recovery by case status (the panel-membership check and the baseline for the drivers);
#   (b) the rare drivers themselves (counts only, in 06_imputed.py);
#   (c) the common-variant leads: the same logistic model on imputed dosages.
# Inputs are read only; no sample IDs are written by 06_imputed.py.
set -euo pipefail
HERE=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck
IMP=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results/04.sample_qc/cteph_agp3k.imputed_array.sqc
MI=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/full_mainland
WGS=$MI/genotype/fixed_model/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready
PHENO=$MI/phenotype/pheno.tsv
COV=$MI/covariates/bbj_mainland_pc.sex.tsv
PL2=/home/b/b37974/plink2
T=${SLURM_CPUS_PER_TASK:-4}
I=$HERE/results/01.inputs
O=$HERE/results/06.imputed; mkdir -p $O

printf '4\t75806990\t76809559\tSTBD1\n12\t112390858\t113390858\tRPH3A\n11\t59076952\t60094182\tOSBP\n9\t4573770\t5573770\tJAK2\n' \
    > $O/drivers.range
printf '%s\n' chr4:76306990:G:C chr4:76309559:G:T chr4:76309487:C:G chr12:112890858:C:T chr11:59576952:G:A \
    chr11:59594182:C:T chr9:5073770:G:T chr9:5073770:G:A > $O/drivers.ids
awk -F'\t' 'NR>1 && $4=="Lead"{print $5}' $I/targets.tsv | sort -u > $O/leads.ids

# imputed: the variants in the driver windows, with their INFO (R2, TYPED/IMPUTED)
awk -F'\t' 'NR==FNR{c[$1]=c[$1]" "$2"-"$3; next} /^#/{next} ($1 in c){n=split(c[$1],w," "); for(i=1;i<=n;i++){split(w[i],b,"-");
    if($2>=b[1] && $2<=b[2]){print $3"\t"$NF; break}}}' $O/drivers.range $IMP.pvar > $O/imp.drivers_window.info
awk -F'\t' 'NR==FNR{a[$1];next} /^#/{next} ($3 in a){print $3"\t"$NF}' $O/leads.ids $IMP.pvar > $O/imp.leads.info
# WGS: allele counts in the 1,683 by case status, in the driver windows
for G in cases controls; do
  $PL2 --bfile $WGS --keep $I/keep.$G.txt --extract range $O/drivers.range --freq counts --out $O/wgs.drv.$G --threads $T
done
# rare in WGS (minor allele count 1-30 in the 1,683) and present in the imputed set; plus the drivers
awk -F'\t' 'FILENAME==ARGV[1]{imp[$1];next}
    FNR==1{for(i=1;i<=NF;i++){if($i=="ALT_CTS")ca=i; if($i=="OBS_CT")co=i}; next}
    FILENAME==ARGV[2]{alt[$2]=$ca; n[$2]=$co; next}
    ($2 in imp){a=alt[$2]+$ca; t=n[$2]+$co; m=(a<t-a)?a:t-a; if(m>=1 && m<=30) print $2}' \
    $O/imp.drivers_window.info $O/wgs.drv.cases.acount $O/wgs.drv.controls.acount > $O/rare_shared.ids
cat $O/rare_shared.ids $O/drivers.ids $O/leads.ids | sort -u > $O/export.ids

$PL2 --bfile $WGS --keep $I/keep.txt --extract $O/export.ids --export A --out $O/wgs --threads $T
$PL2 --pfile $IMP --keep $I/keep.txt --extract $O/export.ids --export A --out $O/imp --threads $T
# the leads on imputed dosages, the same model as steps 1 and 3
$PL2 --pfile $IMP --no-psam-pheno --keep $I/keep.txt --extract $O/leads.ids --pheno $PHENO --pheno-name PHENO1 --1 \
     --covar $COV --covar-name SEX PC1_AVG-PC10_AVG --covar-variance-standardize \
     --glm firth-fallback hide-covar cols=chrom,pos,ref,alt,a1freq,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p \
     --out $O/imp.leads --threads $T
grep -q "11 covariates loaded" $O/imp.leads.log
grep -q "299 cases and 1384 controls remaining" $O/imp.leads.log

/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 -I $HERE/scripts/06_imputed.py
