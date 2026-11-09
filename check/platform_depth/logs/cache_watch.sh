#!/usr/bin/env bash
# Watch a run log for a per-arm task of a genome whose arms were already run (cache expected): any such
# submission is a cache miss; stop the driver by PID and say so. usage: cache_watch.sh <tag>
tag=$1
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/platform_depth
ids=$(paste -sd'|' logs/.cached_ids.txt)
pid=$(cat logs/$tag.pid)
while kill -0 "$pid" 2>/dev/null; do
    hit=$(grep -aE "^\[PROCESS [0-9a-f]{2}/[0-9a-f]{6}\] (EXTRACT_READS|SUBSAMPLE|CALL_GVCF|HLA_SUBSAMPLE|HLAHD) \(($ids)[:)]" logs/$tag.log | head -1)
    if [ -n "$hit" ]; then
        echo "CACHE MISS of an existing arm: $hit; SIGTERM to driver $pid"
        kill "$pid"
        exit 1
    fi
    sleep 5
done
echo "driver $pid exited; no existing arm was re-submitted"
