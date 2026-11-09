#!/usr/bin/env bash
# check/platform_depth — run one stage of platform_depth.nf, detached, behind a preflight and a cache-miss guard.
#
# usage: bash run.sh <stage> [--dry-run] [--queue gr10478b|eb] [-- <more nextflow arguments>]
#        bash run.sh preflight | verify | status [<run tag>]
# stages (arguments and run-tag prefix from scripts/checks.py; docs/PIPELINE.md):
#   inputs  gates  pilot  downsampling  platform  report
#
# preflight (a failure stops the launch):
#   1. no other Nextflow driver of this folder (java ... run platform_depth.nf with this working directory)
#   2. the frozen pieces unchanged: scripts/frozen.py check against docs/FROZEN.manifest.tsv (staged scripts by
#      md5, size and mtime; other staged files; the non-report process definitions; nextflow.config)
#   3. nextflow lint platform_depth.nf nextflow.config (never `run --help`: it launches the pipeline; never -preview)
#   4. the queue: this user's jobs by state and reason (QOSGrpCpuLimit = held by group use: switch with --queue)
# run tag: the next free <prefix>.N (the trace never overwrites, so every run needs a new tag)
# launch: nohup setsid nextflow run ... -resume <session> -ansi-log false > logs/<tag>.log 2>&1 < /dev/null & disown
# guard (detached, logs/<tag>.guard.log): with -ansi-log false a submitted task prints "[PROCESS hh/hhhhhh] NAME" and
#   a cached one prints nothing. A submitted process outside the stage's list (report: the eight report processes;
#   other stages: anything but the inputs processes) is a cache miss: the guard sends SIGTERM to the java driver,
#   by PID (never pkill -f: from a tool shell it matches that shell). GUARD_ALLOW=<regex> overrides the list.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
SESSION=${PD_SESSION:-8eec4a16-3cd8-4479-a6b2-7d6b3d8ecc4c}
PY=/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3
export PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH
export PYTHONDONTWRITEBYTECODE=1

die() { echo "run.sh: $*" >&2; exit 1; }
ts() { date '+%F %T'; }

drivers() {   # PIDs of the Nextflow drivers running platform_depth.nf in this folder [with this run tag]
    local pid args
    for pid in $(ps -u "$(id -u)" -o pid=,comm= | awk '$2 == "java" { print $1 }'); do
        args=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null) || continue
        [[ $args == *" run platform_depth.nf "* ]] || continue
        [[ -z ${1:-} || $args == *"--run_tag $1 "* ]] || continue
        [ "$(readlink "/proc/$pid/cwd" 2>/dev/null)" = "$HERE" ] && echo "$pid"
    done
    return 0
}

next_tag() {  # <prefix> -> <prefix>.N, N above every run of that prefix in results/_run_info and logs
    local p=$1 n=0 f m
    for f in results/_run_info/"$p" results/_run_info/"$p".* logs/"$p".log logs/"$p".*.log; do
        [ -e "$f" ] || continue
        m=${f##*/}; m=${m%.log}; m=${m#"$p"}; m=${m#.}; m=${m%%.*}
        [ -z "$m" ] && m=1
        [[ $m =~ ^[0-9]+$ ]] && (( m > n )) && n=$m
    done
    echo "$p.$((n + 1))"
}

preflight() {
    local d bad=0
    echo "[$(ts)] preflight on $(hostname)"
    d=$(drivers)
    if [ -n "$d" ]; then echo "  FAIL  another driver of this folder is running: PID $d"; bad=1
    else echo "  ok    no other driver of this folder"; fi
    if [ -z "$d" ] && [ -n "$(find .nextflow.log -newermt '2 minutes ago' 2>/dev/null)" ]; then
        echo "  WARN  .nextflow.log changed in the last 2 minutes: a driver on another node?"
    fi
    if "$PY" -B scripts/frozen.py check > logs/.frozen.check 2>&1; then echo "  ok    $(tail -1 logs/.frozen.check)"
    else sed 's/^/  /' logs/.frozen.check; echo "  FAIL  frozen pieces changed (cached tasks would re-run)"; bad=1; fi
    if nextflow lint platform_depth.nf nextflow.config > logs/.lint.out 2>&1; then echo "  ok    nextflow lint"
    else sed 's/^/  /' logs/.lint.out; echo "  FAIL  nextflow lint"; bad=1; fi
    if command -v squeue > /dev/null; then
        q=$(squeue -u "$USER" -h -o '%P:%T:%r' | sort | uniq -c | awk '{ printf "%s %s; ", $1, $2 }')
        echo "  queue ${q:-no jobs of $USER}"
        if squeue -u "$USER" -h -t PD -o '%r' | grep -q QOSGrpCpuLimit; then
            echo "  WARN  jobs held on QOSGrpCpuLimit (group use squeue does not show): consider --queue eb|gr10478b"
        fi
    fi
    return $bad
}

guard() {     # <tag> <driver PID> <allowed process regex>
    local tag=$1 pid=$2 allow=$3 log=logs/$1.log hit
    echo "[$(ts)] guard of $tag: driver $pid; may be submitted: $allow"
    while :; do
        hit=$(grep -aoE '(^\[PROCESS [0-9a-f]{2}/[0-9a-f]{6}\] |Submitted process > )[A-Za-z0-9_:]+' "$log" 2>/dev/null |
              awk '{ print $NF }' | grep -vxE "$allow" | head -1 || true)
        if [ -n "$hit" ]; then
            echo "[$(ts)] CACHE MISS: $hit was submitted; SIGTERM to the driver $pid"
            kill "$pid" 2> /dev/null || true
            for _ in $(seq 60); do kill -0 "$pid" 2> /dev/null || break; sleep 1; done
            kill -0 "$pid" 2> /dev/null && echo "[$(ts)] driver $pid still alive after 60 s: kill it by PID"
            echo "GUARD TRIPPED"
            exit 1
        fi
        kill -0 "$pid" 2> /dev/null || break
        sleep 2
    done
    echo "[$(ts)] driver $pid exited; no cache miss"
    grep -aE '^\[(SUCCESS|ERROR|FAILED)\]|completed=' "$log" | tail -2 || true
    echo "GUARD CLEAN"
}

status() {    # [<tag>]: the last run (or the one named)
    local tag=${1:-$(ls -t logs/*.guard.log 2>/dev/null | head -1 | sed 's#logs/##; s#\.guard\.log$##')}
    [ -n "$tag" ] || die "no guarded run yet"
    echo "run $tag; driver: $(drivers "$tag" | tr '\n' ' ')"
    [ -f "logs/$tag.guard.log" ] && tail -3 "logs/$tag.guard.log"
    [ -f "logs/$tag.log" ] && tail -5 "logs/$tag.log"
    return 0
}

[ $# -ge 1 ] || die "usage: bash run.sh <inputs|gates|pilot|downsampling|platform|report|preflight|verify|status> [...]"
stage=$1; shift
case $stage in
    _guard) guard "$@"; exit $? ;;
    status) status "$@"; exit 0 ;;
    verify) exec bash verify.sh "${1:-downsampling}" ;;
    preflight) preflight; exit $? ;;
esac

spec=$("$PY" -B scripts/checks.py stage "$stage") || die "unknown stage '$stage'"
args=$(sed -n 1p <<< "$spec"); prefix=$(sed -n 2p <<< "$spec"); allow=${GUARD_ALLOW:-$(sed -n 3p <<< "$spec")}
dry=0; extra=()
while [ $# -gt 0 ]; do
    case $1 in
        --dry-run) dry=1 ;;
        --queue) [ $# -ge 2 ] || die "--queue needs a queue"; extra+=(--queue_long "$2" --queue_short "$2"); shift ;;
        --) shift; extra+=("$@"); break ;;
        *) die "unknown option '$1'" ;;
    esac
    shift
done

preflight || die "preflight failed; nothing launched"
tag=$(next_tag "$prefix")
[ ! -e "results/_run_info/$tag" ] && [ ! -e "logs/$tag.log" ] || die "run tag $tag already used"
# shellcheck disable=SC2206
cmd=(nextflow run platform_depth.nf $args --run_tag "$tag" -resume "$SESSION" -ansi-log false "${extra[@]}")
echo "[$(ts)] $tag: ${cmd[*]}"
echo "        guard: may be submitted: $allow"
if [ $dry = 1 ]; then echo "dry run: nothing launched"; exit 0; fi

nohup setsid "${cmd[@]}" > "logs/$tag.log" 2>&1 < /dev/null & disown
pid=''
for _ in $(seq 60); do
    pid=$(drivers "$tag" | head -1)
    [ -n "$pid" ] && break
    sleep 1
done
[ -n "$pid" ] || die "no driver for $tag after 60 s; see logs/$tag.log"
echo "$pid" > "logs/$tag.pid"
nohup setsid bash "$HERE/run.sh" _guard "$tag" "$pid" "$allow" > "logs/$tag.guard.log" 2>&1 < /dev/null & disown
echo "[$(ts)] $tag launched: driver PID $pid (logs/$tag.pid), log logs/$tag.log, guard logs/$tag.guard.log"
echo "        follow with: bash run.sh status $tag"
