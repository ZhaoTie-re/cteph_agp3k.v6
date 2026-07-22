#!/bin/bash
# ---------------------------------------------------------------------------
# Build the association genotypes for ONE replicate x ONE model: the cohort's own
# genotypes with the down-sampled samples' genotypes swapped in.
#
#   $1  REP      replicate label, e.g. rep1
#   $2  MODEL    model label, e.g. snp_based / gene_based
#   $3  BASE     base PLINK prefix (.bed/.bim/.fam) — supplies every sample that
#                was NOT down-sampled, and defines the variant and sample universe
#   $4  REGION   chr:start-end to restrict to
#   $5  POST_QC  the replicate's post-genotype-QC VCF (carries the down-sampled
#                samples only)
#   $6  OUT      output PLINK prefix
#   $7  RELATED  'NONE', or a FID<TAB>IID list of related samples this model must
#                not contain. A fixed-effect model has no GRM to absorb relatedness,
#                so it needs independent samples; a random-effect model keeps them
#                and lets the GRM do the work. Declaring the list turns "the base
#                was already pruned" from an assumption into a checked fact — point
#                a model at the wrong base and it fails here rather than returning
#                a plausible, wrong p-value.
#
# The result is the base fileset with exactly one thing changed: samples that were
# down-sampled carry the replicate's re-called genotypes instead of their original
# ones. Everyone else is untouched. Run the association on this and any effect that
# survives cannot be the depth difference, because it is no longer there.
#
# Variants = REGION ∩ BASE.bim ∩ replicate. The intersection is reported; a base
# variant the replicate never called would otherwise be dropped in silence.
#
# ALLELES. A .bim stores A1/A2, not REF/ALT, so a PLINK->VCF->PLINK round trip can
# quietly swap the two and flip every dosage. The variant IDs are CHROM:POS:REF:ALT,
# so REF is recoverable from the ID itself and is forced at every conversion. The
# script then re-reads the output .bim and fails if any allele disagrees with the
# base — a silent flip here would invert an odds ratio.
#
# Assumes PATH/conda are set by the caller (plink2, bcftools).
# ---------------------------------------------------------------------------
set -euo pipefail

REP=$1; MODEL=$2; BASE=$3; REGION=$4; POST_QC=$5; OUT=$6; RELATED=${7:-NONE}

chr=${REGION%%:*}; range=${REGION#*:}; start=${range%-*}; end=${range#*-}

# ── 0) relatedness: assert, do not assume ─────────────────────────────────
# The base is expected to have been pruned upstream. Expected is not verified, and
# an unpruned base would produce a fixed-effect result that looks entirely normal
# and is quietly wrong, so check before spending anything.
n_related=0
if [ "${RELATED}" != "NONE" ]; then
    awk '{print $2}' "${RELATED}" | sort -u > related_iids.txt
    awk '{print $2}' "${BASE}.fam" | sort -u > base_iids.txt
    comm -12 related_iids.txt base_iids.txt > related_still_in.txt
    n_related=$(wc -l < related_iids.txt)
    n_leak=$(wc -l < related_still_in.txt)
    if [ "${n_leak}" -gt 0 ]; then
        echo "[build_assoc_genotypes] ${REP}/${MODEL}: ${n_leak} related sample(s) are still in" >&2
        echo "  ${BASE}.fam, but this model declares they must not be. Either the base is the" >&2
        echo "  wrong fileset, or it was never pruned:" >&2
        sed 's/^/    /' related_still_in.txt >&2
        exit 1
    fi
    echo "[build_assoc_genotypes] ${REP}/${MODEL}: relatedness checked — none of the ${n_related} related samples are in the base" >&2
fi

# ── 1) variants: the base's own, inside the region ─────────────────────────
# .bim col1 has no 'chr' prefix; col2 is the CHROM:POS:REF:ALT ID; col4 the position.
awk -v c="${chr#chr}" -v s="${start}" -v e="${end}" \
    '$1==c && $4>=s && $4<=e {print $2}' "${BASE}.bim" | sort -u > base_vars.txt

bcftools query -f '%ID\n' "${POST_QC}" | sort -u > rep_vars.txt
comm -12 base_vars.txt rep_vars.txt > vars.txt
comm -23 base_vars.txt rep_vars.txt > vars_not_in_rep.txt

n_base=$(wc -l < base_vars.txt)
n_vars=$(wc -l < vars.txt)
n_gap=$(wc -l < vars_not_in_rep.txt)
if [ "${n_vars}" -eq 0 ]; then
    echo "[build_assoc_genotypes] ${REP}/${MODEL}: no base variant in ${REGION} is in the replicate" >&2
    exit 1
fi
echo "[build_assoc_genotypes] ${REP}/${MODEL}: ${n_base} base variants in ${REGION}, ${n_vars} also in the replicate, ${n_gap} not" >&2

# REF allele from the ID: CHROM:POS:REF:ALT -> field 3.
awk -F: '{print $0"\t"$3}' vars.txt > ref_allele.txt

# ── 2) samples ────────────────────────────────────────────────────────────
# The replicate carries exactly the down-sampled samples, so it defines the swap.
# Emitted alongside the fileset because the downstream comparison needs to know
# which rows were swapped — every other row is identical to the base by
# construction and can only pad a concordance.
bcftools query -l "${POST_QC}" | sort -u > "${OUT}.ds_samples.txt"
awk 'NR==FNR{s[$1];next} ($2 in s){print $1"\t"$2}' "${OUT}.ds_samples.txt" "${BASE}.fam" > ds_keep.txt
n_ds=$(wc -l < ds_keep.txt)
if [ "${n_ds}" -eq 0 ]; then
    echo "[build_assoc_genotypes] ${REP}/${MODEL}: no down-sampled sample is in ${BASE}.fam" >&2
    exit 1
fi
echo "[build_assoc_genotypes] ${REP}/${MODEL}: swapping ${n_ds} of $(wc -l < "${BASE}.fam") samples" >&2

# ── 3) the samples that keep their original genotypes ─────────────────────
# id-paste=iid: the default pastes FID_IID into the VCF sample name, which would
#   not match the replicate's names, and the merge would produce two disjoint
#   sample sets instead of one cohort.
# --output-chr chrM: the .bim numbers chromosomes ('16') while the replicate VCF
#   names them ('chr16'). bcftools merge matches on CHROM:POS:REF:ALT, so the two
#   spellings are two different chromosomes to it: nothing merges, and every
#   variant comes out twice — once per source. Force the VCF spelling here.
plink2 --bfile "${BASE}" \
    --remove ds_keep.txt --extract vars.txt \
    --ref-allele force ref_allele.txt 2 1 \
    --output-chr chrM \
    --export vcf bgz id-paste=iid --out base_nonds
bcftools index -f -t base_nonds.vcf.gz

# ── 4) the down-sampled samples, re-called ────────────────────────────────
bcftools view -i "ID=@vars.txt" "${POST_QC}" -Oz -o rep_part.vcf.gz
bcftools index -f -t rep_part.vcf.gz

# ── 5) one cohort again ───────────────────────────────────────────────────
# The two sample sets are disjoint by construction, so merge cannot collide.
bcftools merge base_nonds.vcf.gz rep_part.vcf.gz -Oz -o merged.vcf.gz
bcftools index -f -t merged.vcf.gz

# --output-chr 26: back to the base's numbering, so the result is a drop-in
# replacement for it rather than something with a differently-spelled chr column.
plink2 --vcf merged.vcf.gz --double-id \
    --ref-allele force ref_allele.txt 2 1 \
    --output-chr 26 \
    --make-bed --out "${OUT}"

n_out=$(wc -l < "${OUT}.bim")
if [ "${n_out}" -ne "${n_vars}" ]; then
    echo "[build_assoc_genotypes] ${REP}/${MODEL}: wrote ${n_out} variants, expected ${n_vars}." >&2
    echo "  A count that is a multiple of the expectation means the merge did not match the" >&2
    echo "  two sources and kept both copies — check CHROM spelling and REF/ALT." >&2
    exit 1
fi

# ── 6) put back what the VCF round trip drops ─────────────────────────────
# A VCF carries no FID, sex or phenotype, so --make-bed invents FID=IID, sex=0,
# pheno=-9. Restore all of it from the base, keyed on IID. Without the phenotype
# the association has nothing to test.
awk 'BEGIN{OFS="\t"}
     NR==FNR {f[$2]=$1; p[$2]=$3; m[$2]=$4; s[$2]=$5; y[$2]=$6; next}
     {
        if (!($2 in f)) { print "sample " $2 " is not in the base fam" > "/dev/stderr"; exit 1 }
        print f[$2], $2, p[$2], m[$2], s[$2], y[$2]
     }' "${BASE}.fam" "${OUT}.fam" > "${OUT}.fam.restored"
mv "${OUT}.fam.restored" "${OUT}.fam"

# ── 7) prove the alleles did not flip ─────────────────────────────────────
# The round trip is the risk; this is the check. Compare A1/A2 against the base
# for every variant that survived.
awk 'BEGIN{OFS="\t"}
     NR==FNR {a1[$2]=$5; a2[$2]=$6; next}
     {
        if (!($2 in a1)) { print "variant " $2 " is not in the base bim" > "/dev/stderr"; bad++; next }
        if ($5 != a1[$2] || $6 != a2[$2]) {
            print "allele mismatch at " $2 ": base " a1[$2] "/" a2[$2] " vs out " $5 "/" $6 > "/dev/stderr"
            bad++
        }
     }
     END { if (bad) { print bad " variant(s) disagree with the base" > "/dev/stderr"; exit 1 } }' \
     "${BASE}.bim" "${OUT}.bim"

# ── 8) summary ────────────────────────────────────────────────────────────
{
    echo "===================================================================================="
    echo " ASSOCIATION GENOTYPES — ${REP} / ${MODEL}"
    echo "===================================================================================="
    echo " generated : $(date -Iseconds)"
    echo " base      : ${BASE}"
    echo " replicate : ${POST_QC}"
    echo " region    : ${REGION}"
    echo " output    : ${OUT}.{bed,bim,fam}"
    echo
    echo " The base cohort with one change: the down-sampled samples carry the replicate's"
    echo " re-called genotypes. Everyone else keeps the genotypes they always had."
    echo
    printf '   %-40s %s\n' "base variants in region" "$(printf "%'d" "${n_base}")"
    printf '   %-40s %s\n' "also called in the replicate -> KEPT" "$(printf "%'d" "${n_vars}")"
    printf '   %-40s %s\n' "not called in the replicate -> dropped" "$(printf "%'d" "${n_gap}")"
    if [ "${n_gap}" -gt 0 ]; then
        echo "     (listed in vars_not_in_rep.txt — these have no re-called genotype to swap in,"
        echo "      so keeping them would mix down-sampled and original data at the same variant)"
    fi
    echo
    printf '   %-40s %s\n' "samples in base" "$(printf "%'d" "$(wc -l < "${BASE}.fam")")"
    printf '   %-40s %s\n' "genotypes swapped for" "$(printf "%'d" "${n_ds}")"
    printf '   %-40s %s\n' "kept as-is" "$(printf "%'d" "$(( $(wc -l < "${BASE}.fam") - n_ds ))")"
    printf '   %-40s %s\n' "samples written" "$(printf "%'d" "$(wc -l < "${OUT}.fam")")"
    printf '   %-40s %s\n' "variants written" "$(printf "%'d" "$(wc -l < "${OUT}.bim")")"
    echo
    echo " CHECKS"
    echo "   alleles vs the base .bim   : all ${n_vars} agree (A1/A2 unflipped)"
    echo "   FID / sex / phenotype      : restored from the base .fam"
    if [ "${RELATED}" != "NONE" ]; then
        echo "   relatedness                : none of the ${n_related} declared related samples are"
        echo "                                present — this model has no GRM to absorb them"
    else
        echo "   relatedness                : not pruned, by design — this model's GRM absorbs it"
    fi
    echo "===================================================================================="
} > "${OUT}.log"

cat "${OUT}.log"
rm -f base_nonds.vcf.gz* rep_part.vcf.gz* merged.vcf.gz*
