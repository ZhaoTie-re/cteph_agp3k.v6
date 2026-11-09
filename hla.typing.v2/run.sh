#!/usr/bin/env bash
# hla.typing.v2 — run every step in order. Plain scripts, no workflow engine: the
# expensive part (HLA-HD over 3,101 samples) was done once by hla.typing and is
# reused, and what remains takes minutes.
#
#   ./run.sh            all steps
#   ./run.sh 01 03      only the named steps
#
# Steps 00a/00b read ~3,100 small archives over Lustre; they run on the login node
# with modest parallelism (12 workers), or submit them with sbatch --rsc if preferred.
set -euo pipefail
cd "$(dirname "$0")"
PY=/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3
mkdir -p logs

steps=("$@")
[ ${#steps[@]} -eq 0 ] && steps=(00a 00b 01 02 03 04 05 06 07)

for s in "${steps[@]}"; do
  case "$s" in
    00a) cmd=(scripts/00_archive_inputs.py --workers 12) ;;
    00b) cmd=(scripts/00_provenance.py) ;;
    01)  cmd=(scripts/01_calls.py) ;;
    02)  cmd=(scripts/02_evidence.py) ;;
    03)  cmd=(scripts/03_panels.py) ;;
    04)  cmd=(scripts/04_confidence.py) ;;
    05)  cmd=(scripts/05_exports.py) ;;
    06)  cmd=(scripts/06_assoc_inputs.py) ;;
    07)  cmd=(scripts/07_figures.py) ;;
    *)   echo "unknown step $s" >&2; exit 2 ;;
  esac
  echo "== step $s: ${cmd[*]}"
  "$PY" -W ignore::FutureWarning "${cmd[@]}" 2>&1 | tee "logs/step_${s}.log"
done
