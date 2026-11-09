#!/usr/bin/env python3
"""Run hla.typing.v2's own steps, unchanged, on another set of typing archives.

hla.typing.v2 keeps every path as a constant in scripts/config.py. This wrapper imports that
module, points RESULTS and every D_* directory at a side tree (ROOT is left alone: it anchors the
pinned IMGT history and fonts), and runs each step script with runpy, so the steps import the
patched module. Nothing is written into hla.typing.v2 (PYTHONDONTWRITEBYTECODE must be set).

  inputs  build a side 00.inputs/: the original archives as symlinks, new archives for the
          re-typed samples (00_archive_inputs.archive_one on the new typing folders), a
          sample_inputs.tsv with their read pairs and archive md5, the pinned freq data
  run     steps 01 -> 02 -> 05 -> 06 (or a subset) on a side tree
  gate    the side tree's 06 matrices against hla.typing.v2's summary.json md5s (gate G3)
"""
import argparse
import gzip
import hashlib
import json
import os
import runpy
import sys
from pathlib import Path

V2 = Path('/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/hla.typing.v2')
STEP_FILES = {'01': '01_calls.py', '02': '02_evidence.py', '05': '05_exports.py', '06': '06_assoc_inputs.py'}


def load_config():
    if os.environ.get('PYTHONDONTWRITEBYTECODE') != '1':
        sys.exit('ABORT: set PYTHONDONTWRITEBYTECODE=1 (no __pycache__ in hla.typing.v2/scripts)')
    sys.path.insert(0, str(V2 / 'scripts'))
    import config as C
    return C


def patch(C, results, inputs):
    results = Path(results)
    C.RESULTS = results
    C.D_INPUTS = Path(inputs)
    C.D_CALLS = results / '01.calls'
    C.D_EVIDENCE = results / '02.evidence'
    C.D_PANELS = results / '03.panels'
    C.D_CONF = results / '04.confidence'
    C.D_EXPORTS = results / '05.exports'
    C.D_ASSOC = results / '06.assoc_inputs'
    C.D_FIG = results / 'figures'
    for d in (C.D_CALLS, C.D_EVIDENCE, C.D_EXPORTS, C.D_ASSOC):
        d.mkdir(parents=True, exist_ok=True)
    return {k: str(getattr(C, k)) for k in ('ROOT', 'RESULTS', 'D_INPUTS', 'D_CALLS', 'D_EVIDENCE',
                                            'D_EXPORTS', 'D_ASSOC')}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fastq_pairs(path):
    n = 0
    with gzip.open(path, 'rt') as fh:
        for _ in fh:
            n += 1
    if n % 4:
        sys.exit(f'ABORT: {path} is not a whole number of FASTQ records')
    return n // 4


def cmd_inputs(a):
    C = load_config()
    orig = V2 / 'results' / '00.inputs'
    side = Path(a.out)
    (side / 'typing').mkdir(parents=True, exist_ok=True)
    retyped = [ln.split('\t')[0] for ln in Path(a.samples).read_text().splitlines() if ln.strip()]
    C.OLD_TYPING = Path(a.typing_dir)
    C.OLD_READS = Path(a.reads_dir)
    sys.argv = ['00_archive_inputs.py']
    arch = runpy.run_path(str(V2 / 'scripts' / '00_archive_inputs.py'), run_name='archive')
    import pandas as pd
    si = pd.read_csv(orig / 'sample_inputs.tsv', sep='\t', dtype=str)
    new = {}
    for sid in retyped:
        _, m, n, _, _ = arch['archive_one'](sid, side / 'typing')
        new[sid] = (m, fastq_pairs(Path(a.reads_dir) / f'{sid}.hla.R1.fastq.gz'))
    for sid in si.sample_id:
        dst = side / 'typing' / f'{sid}.tar.gz'
        if sid not in new and not dst.exists():
            dst.symlink_to(orig / 'typing' / f'{sid}.tar.gz')
    for sid, (m, pairs) in new.items():
        k = si.sample_id == sid
        si.loc[k, 'archive_md5'] = m
        si.loc[k, 'read_pairs'] = str(pairs)
    si.to_csv(side / 'sample_inputs.tsv', sep='\t', index=False)
    fd = side / 'hlahd_freq_data'
    if not fd.exists():
        fd.symlink_to(orig / 'hlahd_freq_data')
    print(f'[hla_v2_wrapper inputs] {len(new)} re-typed archives, {len(si) - len(new)} linked -> {side}')


def cmd_run(a):
    C = load_config()
    used = patch(C, a.results, a.inputs)
    (Path(a.results) / 'wrapper_config.json').write_text(json.dumps(used, indent=2) + '\n')
    for st in a.steps.split(','):
        print(f'[hla_v2_wrapper] step {st} with {used["D_INPUTS"]} -> {used["RESULTS"]}', flush=True)
        sys.argv = [STEP_FILES[st]]
        runpy.run_path(str(V2 / 'scripts' / STEP_FILES[st]), run_name='__main__')


def cmd_gate(a):
    ref = json.loads((V2 / 'results' / '06.assoc_inputs' / 'summary.json').read_text())
    rep, bad = {}, []
    for co, d in ref['cohorts'].items():
        for name, m in d['md5'].items():
            got = md5(Path(a.results) / '06.assoc_inputs' / co / name)
            rep[f'{co}/{name}'] = {'published': m, 'reproduced': got}
            if got != m:
                bad.append(f'{co}/{name}')
    out = {'files': rep, 'mismatches': bad, 'pass': not bad}
    Path(a.out).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'mismatches': bad, 'pass': not bad}))
    if bad:
        sys.exit('GATE G3 FAILED')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('inputs')
    p.add_argument('--samples', required=True, help='the re-typed sample IDs (first column)')
    p.add_argument('--typing-dir', required=True, help='<dir>/<sample>/ = HLA-HD output + run log + validation')
    p.add_argument('--reads-dir', required=True, help='<dir>/<sample>.hla.R1.fastq.gz')
    p.add_argument('--out', required=True)
    p = sub.add_parser('run')
    p.add_argument('--results', required=True)
    p.add_argument('--inputs', required=True)
    p.add_argument('--steps', default='01,02,05,06')
    p = sub.add_parser('gate')
    p.add_argument('--results', required=True)
    p.add_argument('--out', required=True)
    a = ap.parse_args()
    {'inputs': cmd_inputs, 'run': cmd_run, 'gate': cmd_gate}[a.cmd](a)


if __name__ == '__main__':
    main()
