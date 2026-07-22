#!/bin/bash
# ---------------------------------------------------------------------------
# Define the refined_core comparison matrix — the sample x variant set our
# down-sampled genotypes are scored against — once, for every downstream step.
#
#   $1  REF          refined_core PLINK prefix (.bed/.bim/.fam), the truth genotypes
#   $2  CORE_SAMPLES refined_core sample list, "FID <TAB> IID" — defines the cohort
#   $3  OUT_VAR      output: variant IDs, one per line
#   $4  OUT_SMP      output: sample IIDs, one per line
#   $5  OUT_LOG      output: how the matrix was reached
#   $6..             target regions (chr:start-end ...)
#
# SAMPLES   CORE_SAMPLES, as given. This list — not ${REF}.fam — is what defines
#           refined_core, so it is taken at face value and never pre-trimmed here.
#           Each downstream step intersects it with what that step actually has
#           (the replicate VCFs carry only the down-sampled platforms; the PLINK
#           set is relatedness-pruned for the fixed model), and each reports the
#           number it ends up scoring.
# VARIANTS  every ${REF}.bim variant inside the target regions.
#
# Assumes PATH/conda are set by the caller.
# ---------------------------------------------------------------------------
set -euo pipefail

REF=$1; CORE_SAMPLES=$2; OUT_VAR=$3; OUT_SMP=$4; OUT_LOG=$5
shift 5
REGIONS="$*"

# ── variants: target-region rows of the .bim ────────────────────────────────
# .bim col1 is the chromosome WITHOUT the 'chr' prefix (16), while col2, the ID,
# carries it (chr16:53887925:T:C). Match on col1, emit col2.
: > "${OUT_VAR}.tmp"
: > region_counts.txt
for reg in ${REGIONS}; do
    chr=${reg%%:*}; range=${reg#*:}; start=${range%-*}; end=${range#*-}
    n=$(awk -v c="${chr#chr}" -v s="${start}" -v e="${end}" \
        '$1==c && $4>=s && $4<=e {print $2}' "${REF}.bim" | tee -a "${OUT_VAR}.tmp" | wc -l)
    printf '%s\t%s\n' "${reg}" "${n}" >> region_counts.txt
done
sort -u "${OUT_VAR}.tmp" > "${OUT_VAR}"
rm -f "${OUT_VAR}.tmp"

# ── samples: the cohort list, as given ─────────────────────────────────────
cut -f2 "${CORE_SAMPLES}" | sed 's/[[:space:]]*$//' | sort -u > "${OUT_SMP}"

# Reported for context only — never subtracted from the matrix above.
awk '{print $2}' "${REF}.fam" | sort -u > fam_iids.txt
comm -23 "${OUT_SMP}" fam_iids.txt > not_in_plink.txt

n_var=$(wc -l < "${OUT_VAR}")
n_smp=$(wc -l < "${OUT_SMP}")
n_fam=$(wc -l < fam_iids.txt)
n_out=$(wc -l < not_in_plink.txt)

if [ "${n_var}" -eq 0 ]; then
    echo "[build_matrix] no ${REF}.bim variant falls in ${REGIONS}" >&2
    exit 1
fi
if [ "${n_smp}" -eq 0 ]; then
    echo "[build_matrix] ${CORE_SAMPLES} yielded no sample IDs" >&2
    exit 1
fi

{
    echo "===================================================================================="
    echo " REFINED_CORE COMPARISON MATRIX"
    echo "===================================================================================="
    echo " generated : $(date -Iseconds)"
    echo " cohort    : ${CORE_SAMPLES}"
    echo " genotypes : ${REF}"
    echo " regions   : ${REGIONS}"
    echo
    echo " The sample x variant set the replicates are scored on. Both axes are matched"
    echo " BY ID downstream, so anything missing from either side falls out there."
    echo
    echo "------------------------------------------------------------------------------------"
    echo " 1. VARIANTS"
    echo "------------------------------------------------------------------------------------"
    printf '   %-36s %s\n' "region" "variants in .bim"
    while IFS=$'\t' read -r reg n; do
        printf '   %-36s %s\n' "${reg}" "$(printf "%'d" "${n}")"
    done < region_counts.txt
    printf '   %-36s %s\n' "TOTAL (deduplicated)" "$(printf "%'d" "${n_var}")"
    echo
    echo "   From the fixed_model 'fixed_ready' set, which keeps the RARE variants. That is"
    echo "   deliberate: down-sampling removes reads, and the calls most likely to break when"
    echo "   reads are removed are the ones supported by few reads. Scoring only common"
    echo "   variants would hide the effect we are looking for."
    echo
    echo "------------------------------------------------------------------------------------"
    echo " 2. SAMPLES"
    echo "------------------------------------------------------------------------------------"
    printf '   %-36s %s\n' "refined_core (the cohort)" "$(printf "%'d" "${n_smp}")"
    echo
    echo "   Taken from the cohort list as given. refined_core is defined by that list, not"
    echo "   by whichever fileset happens to hold genotypes, so nothing is subtracted here."
    echo "   Each downstream step intersects it with what that step has and reports what it"
    echo "   scored:"
    echo "     - the replicate VCFs carry only the down-sampled platforms"
    printf '     - this PLINK set carries %s of the %s (relatedness-pruned for the\n' \
        "$(printf "%'d" $((n_smp - n_out)))" "$(printf "%'d" "${n_smp}")"
    echo "       fixed model: one of each related pair is dropped so the fixed-effect model"
    echo "       sees independent samples). Those samples are still refined_core."
    if [ "${n_out}" -gt 0 ]; then
        echo
        echo "   refined_core samples with no genotypes in this fileset (${n_out}):"
        sed 's/^/     /' not_in_plink.txt
    fi
    echo "===================================================================================="
} > "${OUT_LOG}"

cat "${OUT_LOG}"
echo "[build_matrix] ${n_var} variants x ${n_smp} refined_core samples" >&2
