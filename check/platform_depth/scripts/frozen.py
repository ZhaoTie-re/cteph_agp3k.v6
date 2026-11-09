#!/usr/bin/env python3
"""The frozen pieces of the cached stages: write or check docs/FROZEN.manifest.tsv (run.sh preflight).

Nextflow's standard cache hashes each task from the process name, its script and every staged file by path +
size + mtime; a change to any of them re-runs the task and everything downstream (~1,000 core-h; a 1 s mtime
change of common.py did once). The manifest records the state the cached runs used:

  script   the scripts staged by cached processes: size, mtime (ns), md5
  staged   other files staged by cached processes (the info sheet, hla.typing.v2's six matrices, the original
           GVCFs JOINT reads): size, mtime; md5 where small
  nf       platform_depth.nf: every process definition but the report stage's, the header (params and helpers;
           comments and blank lines dropped) and the workflow outside its report block: md5 of the text
  config   nextflow.config without comments: md5

  check   exit 1 on any difference, naming it (read-only)
  write   record the current state; only for a deliberate re-freeze (--force to replace a manifest)
"""
import argparse
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs' / 'FROZEN.manifest.tsv'
PROJECT = ROOT.parents[1]
SCRIPTS = ['arm_plan.py', 'bedio.py', 'build_targets.py', 'carrier_firth.py', 'common.py', 'concordance.py',
           'control_subsets.py', 'depth_probe.py', 'depth_summary.py', 'expand_fixed.py', 'gates_ds.py',
           'gene_retest.py', 'genotype_qc.sh', 'hla_concordance.py', 'hla_gate.py', 'hla_matrices.py',
           'hla_v2_wrapper.py', 'make_subsets.py', 'platform_stats.py', 'saige_collect.py', 'sample_table.py',
           'choose_R.py', 't7_fractions.py']
REPORT = ['DS_SUMMARY', 'PLATFORM_REPORT', 'ARRAY_COHORTS', 'ARRAY_EVIDENCE', 'CONTROL_REPORT', 'HLA_REF_SHIFT',
          'FOLLOWUP', 'EVIDENCE', 'FIGURES']       # = checks.REPORT_PROCESSES (verify_checks.py compares them)
SMALL = 200 << 20                      # md5 of staged files up to 200 MB


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as fh:
        for c in iter(lambda: fh.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def text_md5(t):
    return hashlib.md5(t.encode()).hexdigest()


def strip_comments(text):
    """Groovy text without /* */ blocks, // comments outside quotes, and blank lines."""
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
    out = []
    for ln in text.splitlines():
        q, cut = None, len(ln)
        for i, ch in enumerate(ln):
            if q:
                q = None if ch == q else q
            elif ch in '\'"':
                q = ch
            elif ln.startswith('//', i):
                cut = i
                break
        ln = ln[:cut].rstrip()
        if ln.strip():
            out.append(ln)
    return '\n'.join(out) + '\n'


def nf_pieces(path):
    """(item, text) of the frozen parts of the workflow file."""
    lines = Path(path).read_text().splitlines(keepends=True)
    first = next(i for i, ln in enumerate(lines) if ln.startswith('process '))
    pieces = [('platform_depth.nf#header', strip_comments(''.join(lines[:first])))]
    i = 0
    while i < len(lines):
        m = re.match(r'process (\w+) \{', lines[i])
        if m:
            j = next(k for k in range(i + 1, len(lines)) if lines[k].rstrip('\n') == '}')
            if m.group(1) not in REPORT:
                pieces.append((f'platform_depth.nf#{m.group(1)}', ''.join(lines[i:j + 1])))
            i = j + 1
            continue
        if lines[i].startswith('workflow {'):
            body = lines[i:]
            s = next(k for k, ln in enumerate(body) if re.match(r"\s*if \('report' in runs\(\)\)", ln))
            ind = re.match(r'\s*', body[s]).group(0)
            e = next(k for k in range(s + 1, len(body)) if body[k].rstrip('\n') == f'{ind}}}')
            while s > 0 and body[s - 1].strip().startswith('//'):
                s -= 1                         # the report block's own comment lines
            pieces.append(('platform_depth.nf#workflow', strip_comments(''.join(body[:s] + body[e + 1:]))))
            break
        i += 1
    return pieces


def staged_files():
    files = [PROJECT / 'info' / 'cteph_agp3k.v6.20260507.xlsx']
    files += [PROJECT / 'hla.typing.v2' / 'results' / '06.assoc_inputs' / c / f'{k}_dosage.tsv'
              for c in ('narrow_mainland', 'intermediate_mainland', 'full_mainland') for k in ('allele', 'residue')]
    chroms = sorted({ln.split(':')[0] for ln in (ROOT / 'results' / '02.targets' / 'hc.intervals').read_text().split()})
    for label in ('pilot', 'downsampling'):
        f = ROOT / 'results' / '03.calls' / label / 'orig_gvcfs.tsv'
        for ln in f.read_text().splitlines()[1:]:
            sid, d = ln.split('\t')[:2]
            for c in chroms:
                g = Path(d) / f'{sid}.{c}.g.vcf.gz'
                files += [g, Path(f'{g}.tbi')]
    return list(dict.fromkeys(files))


def current():
    rows = []
    for s in SCRIPTS:
        p = ROOT / 'scripts' / s
        st = p.stat()
        rows.append(('script', f'scripts/{s}', str(st.st_size), str(st.st_mtime_ns), md5(p)))
    for p in staged_files():
        if not p.exists():
            rows.append(('staged', str(p), 'missing', '', ''))
            continue
        st = p.stat()
        rows.append(('staged', str(p), str(st.st_size), str(st.st_mtime_ns), md5(p) if st.st_size <= SMALL else ''))
    for item, text in nf_pieces(ROOT / 'platform_depth.nf'):
        rows.append(('nf', item, str(len(text)), '', text_md5(text)))
    cfg = strip_comments((ROOT / 'nextflow.config').read_text())
    rows.append(('config', 'nextflow.config', str(len(cfg)), '', text_md5(cfg)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['check', 'write'])
    ap.add_argument('--manifest', type=Path, default=MANIFEST)
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    rows = current()
    if a.cmd == 'write':
        if a.manifest.exists() and not a.force:
            sys.exit(f'ABORT: {a.manifest} exists; a re-freeze needs --force')
        a.manifest.write_text('kind\titem\tsize\tmtime_ns\tmd5\n' + ''.join('\t'.join(r) + '\n' for r in rows))
        print(f'[frozen] {len(rows)} items -> {a.manifest}')
        return
    want = {ln.split('\t')[1]: ln.split('\t') for ln in a.manifest.read_text().splitlines()[1:]}
    have = {r[1]: list(r) for r in rows}
    bad = []
    for item, w in want.items():
        h = have.get(item)
        if h is None:
            bad.append(f'{item}: gone')
        elif h[2:] != w[2:]:
            what = [k for k, x, y in zip(('size', 'mtime', 'md5'), h[2:], w[2:]) if x != y]
            bad.append(f'{item}: {", ".join(what)} changed')
    bad += [f'{item}: new, not in the manifest' for item in have if item not in want]
    for b in bad:
        print(f'FROZEN CHANGED  {b}')
    print(f'[frozen] {len(want)} items checked, {len(bad)} changed')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
