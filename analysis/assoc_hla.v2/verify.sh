#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# verify.sh — read-only checks of a finished assoc_hla.v2 run.
#
#   §0  the run is independent of analysis/assoc_hla (no path into it)
#   §0  SLURM rules: --rsc only, no cpus/memory directives, m = t x 4571M
#   §0  the archived GEC v0.2 jar is the one the pipeline names (sha256)
#   §1-§13  data checks, recomputed from the inputs (scripts/verify_checks.py)
#
# Exit 0 only if every check passes. Nothing under results/ is written.
# Usage: bash verify.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6
PY=/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3
fail=0

echo '§0 independence from analysis/assoc_hla'
if grep -rnE 'assoc_hla/' "$HERE/assoc_hla.v2.nf" "$HERE/nextflow.config" "$HERE/scripts" >/dev/null 2>&1; then
    echo '  FAIL  a path into analysis/assoc_hla:'
    grep -rnE 'assoc_hla/' "$HERE/assoc_hla.v2.nf" "$HERE/nextflow.config" "$HERE/scripts"
    fail=1
else
    echo '  PASS  no file of this component refers to analysis/assoc_hla'
fi

echo '§0 SLURM resource rules'
if grep -nE '^\s*(cpus|memory)\s' "$HERE/assoc_hla.v2.nf" "$HERE/nextflow.config"; then
    echo '  FAIL  a cpus/memory directive (rejected by the site submitter)'
    fail=1
else
    echo '  PASS  no cpus/memory directive'
fi
# every --rsc request: m must equal t x 4571M (more silently doubles the cores)
bad=0
while read -r t m; do
    [ "$m" -eq $((t * 4571)) ] || { echo "  FAIL  --rsc t=$t m=${m}M is not t x 4571M"; bad=1; }
done < <(grep -oE 't=\$\{params\.[A-Za-z]+\}:c=\$\{params\.[A-Za-z]+\}:m=[0-9]+M' "$HERE/nextflow.config" |
         sed -E 's/t=\$\{params\.([A-Za-z]+)\}.*m=([0-9]+)M/\1 \2/' |
         while read -r p m; do
             t=$(grep -E "^params\.$p\s*=" "$HERE/nextflow.config" | grep -oE '[0-9]+$')
             echo "$t $m"
         done)
[ $bad -eq 0 ] && echo '  PASS  every --rsc memory is cores x 4571M' || fail=1

echo '§0 GEC v0.2 jar'
want=$(grep -oE "GecSha256\s*=\s*'[0-9a-f]{64}'" "$HERE/assoc_hla.v2.nf" | grep -oE '[0-9a-f]{64}')
have=$(sha256sum "$HERE/tools/gec-0.2/gec.jar" | cut -d' ' -f1)
if [ -n "$want" ] && [ "$want" = "$have" ] && grep -q "$have" "$HERE/tools/gec-0.2/SOURCE.md"; then
    echo "  PASS  tools/gec-0.2/gec.jar sha256 = params.GecSha256 = SOURCE.md ($have)"
else
    echo "  FAIL  gec.jar sha256 $have; params.GecSha256 '$want'"
    fail=1
fi

"$PY" "$HERE/scripts/verify_checks.py" \
    --results "$HERE/results" \
    --typing "$PROJECT/hla.typing.v2/results/06.assoc_inputs" \
    --model-inputs "$PROJECT/wgs.auto.par/results/12_model_inputs" \
    --gwas-dir "$PROJECT/analysis/assoc_saige/results" \
    --gec-jar "$HERE/tools/gec-0.2/gec.jar" \
    --tmp "$(mktemp -d "${TMPDIR:-/tmp}/assoc_hla_v2_verify.XXXXXX")" || fail=1

if [ $fail -eq 0 ]; then echo 'verify.sh: ALL CHECKS PASSED'; else echo 'verify.sh: FAILURES ABOVE'; fi
exit $fail
