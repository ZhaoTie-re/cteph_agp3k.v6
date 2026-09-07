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
#           --crai-dir searches a third place last: a directory of indexes this
#           pipeline built earlier. IT IS NOT WIRED UP, deliberately — see the note
#           on BUILD_MANIFEST in hla.typing.nf. Pointing it at the run's own output
#           makes the manifest a function of the previous run rather than of the
#           inputs, and the second run then re-extracts and re-types the 111 samples
#           whose crai column changed. The option stays for a caller who has an
#           index directory that is genuinely an INPUT.
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
    p.add_argument('--cohort-keep', default=None,
                   help='the analysis cohort, e.g. PopGMM_output/full_mainland.fid_iid.txt. '
                        'Applied AFTER --keep-id and kept separate from it on purpose: '
                        'sample QC decides what is analysable, the cohort list decides what '
                        'is COMPARABLE, and collapsing the two would lose the distinction. '
                        'Same FID<TAB>IID shape. Omit to type every QC-passing sample.')
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


def read_keep(path, what):
    """The ids in a FID<TAB>IID keep list, or a named abort.

    Both lists this component reads — wgs.auto.par's sample_qc.keep.id and
    PopGMM_output/<cohort>.fid_iid.txt — are the same two-column shape with the
    sample id in both columns, so one reader serves both.
    """
    path = Path(path)
    if not path.exists():
        raise SystemExit(
            f'ABORT: no {what} at {path}.\n'
            f'  The upstream component has not produced it, or its work/ was cleaned.')
    ids = {ln.split()[0] for ln in path.read_text().splitlines() if ln.strip()}
    if not ids:
        raise SystemExit(f'ABORT: {path} is empty')
    return ids


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
    keep = read_keep(keep_path, 'sample-QC keep list (wgs.auto.par 07_sample_qc)')

    sel = sel[sel[args.id_col].astype(str).isin(keep)]
    n_qc_dropped = n_flag - len(sel)
    if not len(sel):
        raise SystemExit(
            f'ABORT: none of the {n_flag} flagged samples is in {keep_path}. '
            f'The two are probably keyed on different ids.')
    n_qc_pass = len(sel)

    # THE THIRD STAGE, and it is an ANCESTRY filter, not a quality one. Sample QC
    # leaves 3,569 analysable samples; the cohort list is the subset that shares an
    # ancestry background, so that an allele absent from a Japanese reference panel
    # means something about the typing rather than about the sample. Typing every
    # QC-passing sample and subsetting later is the other valid choice; this
    # component takes the cohort as its scope so that every table it publishes has
    # one denominator.
    n_cohort_dropped = 0
    cohort_path = Path(args.cohort_keep) if args.cohort_keep else None
    if cohort_path is not None:
        cohort = read_keep(cohort_path, 'cohort keep list')
        absent = sorted(cohort - set(sel[args.id_col].astype(str)))
        if absent:
            raise SystemExit(
                f'ABORT: {len(absent)} of {len(cohort)} cohort sample(s) in {cohort_path} '
                f'are not available to type — they are absent from {args.xlsx} with '
                f'{args.flag_col} == True, or were dropped by {keep_path.name}.\n'
                f'  first few: {absent[:10]}\n'
                f'  The cohort must be a SUBSET of the QC-passing set; if it is not, the '
                f'two were built from different sample universes.')
        sel = sel[sel[args.id_col].astype(str).isin(cohort)]
        n_cohort_dropped = n_qc_pass - len(sel)

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
    print(f'    {args.flag_col} == True: {n_flag}')
    print(f'    -{n_qc_dropped:<5d} sample QC ({keep_path.name})  -> {n_qc_pass}')
    if cohort_path is not None:
        print(f'    -{n_cohort_dropped:<5d} cohort   ({cohort_path.name})  -> {len(sel)}')
    else:
        print('    (no --cohort-keep: every QC-passing sample is typed)')
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
