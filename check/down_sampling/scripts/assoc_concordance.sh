#!/bin/bash
# ---------------------------------------------------------------------------
# Score every association fileset against its base, on the variants that test
# actually reads. Wraps assoc_concordance.py.
#
#   $1  PY       path to assoc_concordance.py
#   $2  SPEC     one line per comparison: MODEL <TAB> REP <TAB> BASE_PREFIX
#                (the assoc filesets are expected in the CWD as <REP>.<MODEL>.bed)
#
# Both sides are exported as ALT dosage with the ALT named explicitly, so the
# comparison cannot be flipped by whichever allele PLINK would otherwise count.
# Assumes PATH/conda are set by the caller (plink2, python3).
# ---------------------------------------------------------------------------
set -euo pipefail

PY=$1; SPEC=$2

: > ds_samples.txt
PAIRS=()
while IFS=$'\t' read -r model rep base; do
    [ -z "${model:-}" ] && continue
    pfx="${rep}.${model}"

    # The variants this model reads, and the ALT to count for each. The ID is
    # CHROM:POS:REF:ALT, so the ALT is the last colon-separated field — naming it
    # keeps the dosage ALT-relative on both sides.
    awk '{print $2}' "${pfx}.bim" | sort -u > "${pfx}.vars.txt"
    awk -F: '{print $0"\t"$NF}' "${pfx}.vars.txt" > "${pfx}.alt.txt"

    # Same samples, same variants, same counted allele — only the genotypes differ.
    plink2 --bfile "${base}" --keep "${pfx}.fam" --extract "${pfx}.vars.txt" \
        --export A --export-allele "${pfx}.alt.txt" --out "${pfx}.base"
    plink2 --bfile "${pfx}" --extract "${pfx}.vars.txt" \
        --export A --export-allele "${pfx}.alt.txt" --out "${pfx}.assoc"

    PAIRS+=(--pair "${model}:${rep}:${pfx}.base.raw:${pfx}.assoc.raw")
done < "${SPEC}"

if [ ${#PAIRS[@]} -eq 0 ]; then
    echo "[assoc_concordance] ${SPEC} listed no comparison" >&2
    exit 1
fi

# The down-sampled samples: whoever the replicate VCFs carried. Every other sample
# is byte-identical between base and assoc, so scoring them would only pad the
# diagonal with rows that cannot disagree.
cat *.ds_samples.txt | sort -u > ds_samples.txt
if [ ! -s ds_samples.txt ]; then
    echo "[assoc_concordance] no down-sampled sample list found" >&2
    exit 1
fi
echo "[assoc_concordance] scoring $(wc -l < ds_samples.txt) down-sampled samples" >&2

python3 "${PY}" \
    "${PAIRS[@]}" \
    --ds-samples ds_samples.txt \
    --out-log assoc_concordance.log \
    --out-fig assoc_concordance.png
