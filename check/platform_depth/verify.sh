#!/usr/bin/env bash
# Read-only verification of check/platform_depth (exit 0 = every check passed).
# usage: bash verify.sh [run label, default downsampling]
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3 scripts/verify_checks.py --label "${1:-downsampling}"
