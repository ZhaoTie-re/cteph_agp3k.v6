#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_calc_metrics.sh [options]

Required:
  --bed-prefix PATH
  --info-file PATH
  --out-sample FILE
  --out-variant FILE
  --script-dir PATH
  --plink2 PATH
  --tabix PATH
  --id-col NAME
  --group-col NAME
  --case-value VALUE
  --tdp-col NAME
  --mdp-col NAME

Optional:
  --min-ac INT      Minor allele count cutoff (default: 0)
  --threads INT     Number of threads (default: 1)
EOF
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

log() {
    echo "[$(date '+%F %T')] $*"
}

require_file_prefix() {
    local prefix="$1"
    [[ -f "${prefix}.bed" && -f "${prefix}.bim" && -f "${prefix}.fam" ]] || {
        die "PLINK prefix not complete: ${prefix}.{bed,bim,fam}"
    }
}

THREADS=1
MIN_AC=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --bed-prefix) BED_PREFIX="${2:-}"; shift 2 ;;
        --info-file) INFO_FILE="${2:-}"; shift 2 ;;
        --out-sample) OUT_SAMPLE="${2:-}"; shift 2 ;;
        --out-variant) OUT_VARIANT="${2:-}"; shift 2 ;;
        --script-dir) SCRIPT_DIR="${2:-}"; shift 2 ;;
        --plink2) PLINK2="${2:-}"; shift 2 ;;
        --tabix) TABIX="${2:-}"; shift 2 ;;
        --id-col) ID_COL="${2:-}"; shift 2 ;;
        --group-col) GROUP_COL="${2:-}"; shift 2 ;;
        --case-value) CASE_VALUE="${2:-}"; shift 2 ;;
        --tdp-col) TDP_COL="${2:-}"; shift 2 ;;
        --mdp-col) MDP_COL="${2:-}"; shift 2 ;;
        --min-ac) MIN_AC="${2:-}"; shift 2 ;;
        --threads) THREADS="${2:-}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown parameter: $1" ;;
    esac
done

[[ -n "${BED_PREFIX:-}" ]] || { usage; die "--bed-prefix is required"; }
[[ -n "${INFO_FILE:-}" ]] || { usage; die "--info-file is required"; }
[[ -n "${OUT_SAMPLE:-}" ]] || { usage; die "--out-sample is required"; }
[[ -n "${OUT_VARIANT:-}" ]] || { usage; die "--out-variant is required"; }
[[ -n "${SCRIPT_DIR:-}" ]] || { usage; die "--script-dir is required"; }
[[ -n "${PLINK2:-}" ]] || { usage; die "--plink2 is required"; }
[[ -n "${TABIX:-}" ]] || { usage; die "--tabix is required"; }
[[ -n "${ID_COL:-}" ]] || { usage; die "--id-col is required"; }
[[ -n "${GROUP_COL:-}" ]] || { usage; die "--group-col is required"; }
[[ -n "${CASE_VALUE:-}" ]] || { usage; die "--case-value is required"; }
[[ -n "${TDP_COL:-}" ]] || { usage; die "--tdp-col is required"; }
[[ -n "${MDP_COL:-}" ]] || { usage; die "--mdp-col is required"; }

[[ "$MIN_AC" =~ ^[0-9]+$ ]] || die "--min-ac must be a non-negative integer"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"

[[ -x "$PLINK2" ]] || die "plink2 not executable: $PLINK2"
[[ -x "$TABIX" ]] || die "tabix not executable: $TABIX"
PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="python3"
    elif command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
    else
        die "Neither python3 nor python is available in PATH"
    fi
fi

command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python executable not found: $PYTHON_BIN"
command -v bgzip >/dev/null 2>&1 || die "bgzip not found in PATH"
[[ -f "$SCRIPT_DIR/calc_metrics.py" ]] || die "Missing script: $SCRIPT_DIR/calc_metrics.py"
[[ -f "$INFO_FILE" ]] || die "Info file not found: $INFO_FILE"
require_file_prefix "$BED_PREFIX"

tmp_dir="$(mktemp -d -p "$PWD" calc_metrics.XXXXXX)"
cleanup() {
    rm -rf "$tmp_dir"
}
trap cleanup EXIT

log "Starting metrics calculation"
log "BED prefix: $BED_PREFIX"
log "MinAC cutoff: $MIN_AC"
log "Threads: $THREADS"

current_bed="$BED_PREFIX"

if [[ "$MIN_AC" -gt 0 ]]; then
    log "Step 0/6: Applying MAC filter >= $MIN_AC"
    filtered_prefix="$tmp_dir/filtered.minac${MIN_AC}"
    "$PLINK2" --threads "$THREADS" --bfile "$BED_PREFIX" --mac "$MIN_AC" --make-bed --out "$filtered_prefix" --silent
    require_file_prefix "$filtered_prefix"
    current_bed="$filtered_prefix"
else
    log "Step 0/6: Skip filtering because MinAC = 0"
fi

log "Step 1/6: Sample counts and missingness"
"$PLINK2" --threads "$THREADS" --bfile "$current_bed" --sample-counts cols=maybefid,homref,het,homalt --missing sample-only --out "$tmp_dir/temp_sample"

log "Step 2/6: Sample minor allele burden (SMinAC)"
"$PLINK2" --threads "$THREADS" --bfile "$current_bed" --maj-ref force --make-bed --out "$tmp_dir/temp_aligned"
"$PLINK2" --threads "$THREADS" --bfile "$tmp_dir/temp_aligned" --sample-counts cols=maybefid,het,homalt --out "$tmp_dir/temp_sample_minac"

log "Step 3/6: Variant counts and missingness"
"$PLINK2" --threads "$THREADS" --bfile "$current_bed" --geno-counts cols=chrom,pos,ref,alt,homref,refalt1,homalt1 --missing variant-only --out "$tmp_dir/temp_variant"

log "Step 4/6: Allele count table"
"$PLINK2" --threads "$THREADS" --bfile "$current_bed" --freq counts --out "$tmp_dir/temp_freq"

log "Step 5/6: Aggregate metrics with Python"
"$PYTHON_BIN" "$SCRIPT_DIR/calc_metrics.py" \
    --sample-counts "$tmp_dir/temp_sample.scount" \
    --sample-missing "$tmp_dir/temp_sample.smiss" \
    --sample-minac "$tmp_dir/temp_sample_minac.scount" \
    --variant-counts "$tmp_dir/temp_variant.gcount" \
    --variant-missing "$tmp_dir/temp_variant.vmiss" \
    --freq-counts "$tmp_dir/temp_freq.acount" \
    --info "$INFO_FILE" \
    --id-col "$ID_COL" \
    --group-col "$GROUP_COL" \
    --case-value "$CASE_VALUE" \
    --tdp-col "$TDP_COL" \
    --mdp-col "$MDP_COL" \
    --out-sample "$OUT_SAMPLE" \
    --out-variant "$OUT_VARIANT" \
    --threads "$THREADS" \
    --log calc_metrics.log

log "Step 6/6: bgzip + tabix outputs"
bgzip -f "$OUT_SAMPLE"
bgzip -f "$OUT_VARIANT"
"$TABIX" -f -s 1 -b 2 -e 2 "${OUT_VARIANT}.gz"

log "Completed successfully"
