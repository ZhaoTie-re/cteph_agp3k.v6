#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Provenance of one snapshot, written once its workbook exists:
#             run_manifest.json  the run's parameters and script md5s (from the .nf),
#                                the Open Targets release, the record counts, the
#                                workbook's rows per sheet and the change summary
#             latest.json        the newest complete snapshot under the results root
#                                (by fetch date), for the knowledge-base pointer
# Inputs  : the run json from the .nf, tidy.json, workbook.json, changes.json
# Outputs : run_manifest.json, latest.json
# ---------------------------------------------------------------------------
import argparse
import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ot_common import csv_list  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-json', required=True, help='parameters and script md5s, from the .nf')
    ap.add_argument('--tidy', required=True)
    ap.add_argument('--workbook-stats', required=True)
    ap.add_argument('--changes-json', required=True)
    ap.add_argument('--snapshot', required=True, help="this snapshot's name (its fetch date)")
    ap.add_argument('--earlier', default='', help='the other snapshots that exist (comma-separated)')
    ap.add_argument('--snapshots-root', required=True)
    ap.add_argument('--workbook', required=True, help="the workbook's file name")
    ap.add_argument('--out', default='run_manifest.json')
    ap.add_argument('--out-latest', default='latest.json')
    a = ap.parse_args()
    run = json.loads(Path(a.run_json).read_text())
    tidy = json.loads(Path(a.tidy).read_text())
    wb = json.loads(Path(a.workbook_stats).read_text())
    ch = json.loads(Path(a.changes_json).read_text())
    man = {'component': 'cteph_knowledge_graph', 'snapshot': a.snapshot, 'written_at':
           dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'), 'open_targets': {
               'data_version': tidy['release'], 'api_version': tidy['api_version'], 'fetched_at': tidy['fetched_at']},
           'counts': {k: v for k, v in tidy.items() if k not in ('release', 'api_version', 'fetched_at')},
           'workbook': {'file': a.workbook, 'rows': wb['rows']}, 'changes': ch, **run}
    Path(a.out).write_text(json.dumps(man, indent=1))
    names = sorted(set(csv_list(a.earlier)) | {a.snapshot})
    newest = names[-1]
    root = Path(a.snapshots_root)
    latest = {'snapshot': newest, 'dir': str(root / newest),
              'knowledge': str(root / newest / '06.knowledge' / 'cteph_knowledge_graph.md'),
              'snapshots': names}
    if newest == a.snapshot:
        latest.update({'release': tidy['release'], 'workbook': str(root / newest / a.workbook)})
    else:   # an older date was re-run: the pointer stays on the newer snapshot, as its own manifest records it
        other = json.loads((root / newest / '_run_info' / 'run_manifest.json').read_text())
        latest.update({'release': other['open_targets']['data_version'],
                       'workbook': str(root / newest / other['workbook']['file'])})
    Path(a.out_latest).write_text(json.dumps(latest, indent=1))


if __name__ == '__main__':
    main()
