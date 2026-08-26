#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Turn the master sample sheet into a per-sample manifest of CRAMs and
#           their indexes, and refuse to produce one at all if any sample is not
#           usable.
#
#           The index is resolved through realpath, not by appending '.crai' to
#           the path in the sheet. 3,117 of this study's 3,569 CRAM paths are
#           symlinks whose index sits beside the REAL file, so the appended form
#           does not exist and samtools fails with
#           `[E::cram_index_load] Could not retrieve index file`.
#
#           A third place is searched last: the run's own crai/ directory. 111 CRAMs
#           have no index anywhere, INDEX_CRAM builds one for each — ~30 minutes of
#           CRAM reading apiece — and publishes it there. Without looking, a run that
#           lost work/ rebuilds all 111 from scratch.
#
#           Failing here is the point. A sample dropped silently at this stage is
#           a sample missing from the cohort with nothing in the logs to say so.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# Index spellings in use, most specific first. samtools writes <f>.cram.crai by
# default but <f>.crai is equally valid and both appear in shared trees.
INDEX_SUFFIXES = ('.crai', '.cram.crai')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xlsx', required=True, help='master sample sheet')
    p.add_argument('--id-col', default='ID_JHRPv6')
    p.add_argument('--flag-col', default='Flag_JHRPv6')
    p.add_argument('--cram-col', default='Cram_Path')
    p.add_argument('--platform-col', default='WGS_Platform')
    p.add_argument('--group-col', default='Outcome')
    p.add_argument('--sex-col', default='Sex')
    p.add_argument('--depth-col', default='Observed_Depth')
    p.add_argument('--keep-id', required=True,
                   help='sample_qc.keep.id from wgs.auto.par — the samples that '
                        'passed sample QC. FID<TAB>IID, both the sample id.')
    p.add_argument('--crai-dir', default=None,
                   help="this run's own crai/, searched after the two beside-the-CRAM "
                        'spellings. Indexes INDEX_CRAM built on an earlier run are '
                        'reused instead of rebuilt.')
    p.add_argument('--max-samples', type=int, default=0,
                   help='0 = all. A positive value takes the first N of each '
                        'platform in turn, so a pilot spans every platform '
                        'rather than the first N rows of one.')
    p.add_argument('--out', default='sample_manifest.tsv')
    return p.parse_args()


def resolve_index(cram, sample=None, crai_dir=None):
    """The index for `cram`, following the symlink to where it actually lives.

    Beside the real file first — both spellings — then in this run's own crai/, where
    INDEX_CRAM publishes the ones it had to build. Returns None when there is none,
    so the caller can name the sample.
    """
    real = os.path.realpath(cram)
    stem = real[:-5] if real.endswith('.cram') else real
    cands = [real + '.crai', stem + '.crai']
    if crai_dir and sample:
        cands.append(os.path.join(crai_dir, f'{sample}.cram.crai'))
    for cand in cands:
        if os.path.exists(cand):
            return cand
    return None


def main():
    args = parse_args()
    d = pd.read_excel(args.xlsx, engine='openpyxl', dtype={args.id_col: str})

    # Every column the manifest carries, not only the four the workflow cannot run
    # without. `group`, `sex` and `depth` were read with .get() and became silently
    # empty when absent, which turns a renamed column into a QC table with no rows
    # rather than into an error.
    for col in (args.id_col, args.flag_col, args.cram_col, args.platform_col,
                args.group_col, args.sex_col, args.depth_col):
        if col not in d.columns:
            raise SystemExit(f'ABORT: {args.xlsx} has no column {col!r}. '
                             f'Found: {list(d.columns)}')

    sel = d[d[args.flag_col] == True].copy()          # noqa: E712 — a real bool column
    if not len(sel):
        raise SystemExit(f'ABORT: no row has {args.flag_col} == True in {args.xlsx}')
    n_flag = len(sel)

    # The workbook flag defines the v6 cohort; SAMPLE QC decides which of it is
    # analysable, and this component inherits that decision rather than making its
    # own. The exclusions are heterozygosity outliers, which for HLA typing are
    # exactly the samples least likely to resolve into two clean haplotypes.
    keep_path = Path(args.keep_id)
    if not keep_path.exists():
        raise SystemExit(
            f'ABORT: no sample-QC keep list at {keep_path}.\n'
            f'  wgs.auto.par has not produced 07_sample_qc, or its work/ was cleaned.\n'
            f'  Run wgs.auto.par first, or point --keep-id at the list to use.')
    keep = {ln.split()[0] for ln in keep_path.read_text().splitlines() if ln.strip()}
    if not keep:
        raise SystemExit(f'ABORT: {keep_path} is empty')

    sel = sel[sel[args.id_col].astype(str).isin(keep)]
    n_qc_dropped = n_flag - len(sel)
    if not len(sel):
        raise SystemExit(
            f'ABORT: none of the {n_flag} flagged samples is in {keep_path}. '
            f'The two are probably keyed on different ids.')

    # Every problem is collected before anything is reported, so one run names all
    # of them instead of one per re-launch.
    no_id = sel[sel[args.id_col].isna()][args.id_col].index.tolist()
    no_cram = sel[sel[args.cram_col].isna()]
    if len(no_id) or len(no_cram):
        raise SystemExit(
            f'ABORT: {len(no_id)} row(s) with no {args.id_col} and '
            f'{len(no_cram)} with no {args.cram_col} among {len(sel)} selected.\n'
            f'  first few missing a CRAM: {list(no_cram[args.id_col].head(10))}')

    rows, missing_cram = [], []
    for _, r in sel.iterrows():
        sid, cram = str(r[args.id_col]), str(r[args.cram_col])
        if not os.path.exists(cram):
            missing_cram.append((sid, cram))
            continue
        # An empty `crai` is not a failure: the workflow indexes those samples into
        # its own tree first. A missing CRAM is a failure — nothing can recover it.
        rows.append({
            'sample_id': sid,
            'cram': cram,
            'crai': resolve_index(cram, sid, args.crai_dir) or '',
            'platform': r[args.platform_col],
            'group': r.get(args.group_col, ''),
            'sex': r.get(args.sex_col, ''),
            'observed_depth': r.get(args.depth_col, ''),
        })

    if missing_cram:
        msg = [f'ABORT: {len(missing_cram)} unreadable CRAM(s) of {len(sel)} selected.']
        for sid, path in missing_cram[:10]:
            msg.append(f'  {sid}  {path}')
        raise SystemExit('\n'.join(msg))

    out = pd.DataFrame(rows)
    if args.max_samples > 0:
        # Round-robin across platforms so a pilot exercises every one of them.
        # Round-robin across platforms, LARGEST FIRST, and interleaved rather than
        # concatenated. Taking each platform's first k and then truncating drops the
        # last platforms alphabetically — which for this study meant a 2-sample pilot
        # was two cases and no control, so it could not exercise anything that needs
        # a control group. Largest first also means the platform carrying most of the
        # cohort is the one a small pilot is guaranteed to test.
        order = out.platform.value_counts().index.tolist()
        rank = {pf: i for i, pf in enumerate(order)}
        out = out.assign(_p=out.platform.map(rank),
                         _i=out.groupby('platform').cumcount())
        out = (out.sort_values(['_i', '_p'])
                  .head(args.max_samples)
                  .drop(columns=['_p', '_i']))
    out = out.sort_values('sample_id').reset_index(drop=True)
    out.to_csv(args.out, sep='\t', index=False)

    n_noidx = int((out.crai == '').sum())
    print(f'[build_manifest] {len(out)} sample(s) -> {args.out}')
    print(f'    {args.flag_col} == True: {n_flag}; '
          f'{n_qc_dropped} dropped by sample QC ({keep_path.name}); {n_flag - n_qc_dropped} kept')
    for pf, g in out.groupby('platform'):
        print(f'    {str(pf):22s} {len(g):5d}')
    print(f'    symlinked CRAMs (index resolved via realpath): '
          f'{sum(os.path.islink(c) for c in out.cram)}')
    if args.crai_dir:
        n_reused = int(sum(str(c).startswith(str(args.crai_dir)) for c in out.crai))
        print(f'    indexes reused from {args.crai_dir}: {n_reused}')
    print(f'    CRAMs with NO index anywhere, to be indexed by this run: {n_noidx}')
    if n_noidx:
        for _, r in out[out.crai == ''].head(5).iterrows():
            print(f'        {r.sample_id}  {os.path.realpath(r.cram)}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in build_manifest: {e}', file=sys.stderr)
        sys.exit(1)
