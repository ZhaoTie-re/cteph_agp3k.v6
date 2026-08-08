#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Purpose : Shared structured logging for tuning.rv Nextflow processes. Every
#           stage sources this and gets one consistent, human-readable record:
#           a banner header, `[STAGE | ISO-8601] key : value` lines, all
#           written to BOTH the console (.command.out) and the published log.
# Project : cteph_agp3k.v6 / tuning.rv  (rare-variant depth-confounding QC)
# Usage   : source rv_log.sh ; rv_log_init CALC_METRICS calc_metrics.log
#           rv_kv "MinAC" 5 ; rv_log "step 1/8 ..." ; rv_done
# ---------------------------------------------------------------------------
RV_STAGE="${RV_STAGE:-STAGE}"
RV_LOG="${RV_LOG:-/dev/null}"

_rv_ts() { date '+%Y-%m-%dT%H:%M:%S%z'; }

# rv_log_init <STAGE> <LOGFILE> — start a fresh log with a banner header.
rv_log_init() {
    RV_STAGE="$1"
    RV_LOG="$2"
    : > "$RV_LOG"
    {
        echo "==================================================================="
        echo "  tuning.rv · ${RV_STAGE}"
        echo "  started : $(_rv_ts)"
        echo "  host    : $(hostname 2>/dev/null || echo '?')"
        echo "  workdir : $(pwd)"
        echo "==================================================================="
    } | tee -a "$RV_LOG"
}

# rv_log <message> — timestamped, stage-tagged line to console + log.
rv_log() { echo "[${RV_STAGE} | $(_rv_ts)] $*" | tee -a "$RV_LOG"; }

# rv_kv <key> <value> — aligned key/value provenance line.
rv_kv() { printf '[%s | %s] %-24s : %s\n' "$RV_STAGE" "$(_rv_ts)" "$1" "$2" | tee -a "$RV_LOG"; }

# rv_section <title> — a labelled divider within the log.
rv_section() {
    { echo ""; echo "--- $* ---"; } | tee -a "$RV_LOG"
}

# rv_run <cmd...> — echo then run a command, appending its output to the log only.
rv_run() {
    rv_log "\$ $*"
    "$@" >> "$RV_LOG" 2>&1
}

# rv_done — closing line.
rv_done() { rv_log "completed successfully."; }
