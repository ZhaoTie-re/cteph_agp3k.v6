#!/usr/bin/env python3
"""Step 00b — what was actually run, read back from the archived task records.

For every sample this records, from its archive:
  - the HLA-HD version printed in its own run log;
  - the HLA-HD command (sample name masked) and the read-extraction command;
  - read pairs written for HLA-HD and singletons discarded (extraction log);
  - the validation status written when the typing ran.

The run is accepted as one consistent typing only if every sample shows the same
version and the same two commands and passed validation; otherwise this step stops.

Writes results/00.inputs/sample_inputs.tsv (one row per sample, with group, sex,
platform and depth from the earlier manifest) and results/00.inputs/provenance.json
(versions, commands, md5 of every pinned file, counts).
"""
import hashlib
import json
import re
import shutil
import sys
import tarfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

VER_RE = re.compile(r'HLA-HD version ([0-9.]+)')
PAIRS_RE = re.compile(r':\s*([\d,]+) read pairs over (\d+) region')
SINGLE_RE = re.compile(r'discarded (\d+) singletons')
PROC_RE = re.compile(r'processed (\d+) reads')


def md5_file(p):
    h = hashlib.md5()
    with open(p, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def md5_dir(d):
    """One digest over every file's relative path and md5, sorted."""
    h = hashlib.md5()
    for p in sorted(x for x in Path(d).rglob('*') if x.is_file()):
        h.update(f'{p.relative_to(d)}\t{md5_file(p)}\n'.encode())
    return h.hexdigest()


def mask(cmd, sid):
    lines = [ln.strip() for ln in cmd.splitlines()
             if ln.strip() and not ln.strip().startswith('#')]
    return '\n'.join(lines).replace(sid, '<SAMPLE>')


def one(path):
    sid = path.name[:-len('.tar.gz')]
    files = {}
    with tarfile.open(path, 'r:gz') as tf:
        for m in tf.getmembers():
            n = m.name.split('/', 1)[1]
            if n in ('hlahd.run.log', 'typing_validation.tsv') or n.startswith('task/'):
                files[n] = tf.extractfile(m).read().decode('utf-8', errors='replace')
    ver = VER_RE.search(files.get('hlahd.run.log', ''))
    out = files.get('task/extract.command.out', '')
    err = files.get('task/extract.command.err', '')
    mp, ms, mr = PAIRS_RE.search(out), SINGLE_RE.search(err), PROC_RE.search(err)
    val = files.get('typing_validation.tsv', '').strip().splitlines()
    status = val[1].split('\t')[-1] if len(val) > 1 else ''
    # The HLA-HD call itself: comments dropped first, then the lines from the
    # `.../hlahd.sh` invocation up to its output redirection.
    hl = mask(files.get('task/hlahd.command.sh', ''), sid)
    a = hl.find('/bin/hlahd.sh')
    a = hl.rfind('\n', 0, a) + 1 if a >= 0 else 0
    b = hl.find('> hlahd.run.log')
    hl_cmd = hl[a:b if b > a else len(hl)].strip()
    return {
        'sample_id': sid,
        'hlahd_version': ver.group(1) if ver else '',
        'read_pairs': int(mp.group(1).replace(',', '')) if mp else -1,
        'regions': int(mp.group(2)) if mp else -1,
        'singletons_discarded': int(ms.group(1)) if ms else -1,
        'reads_processed': int(mr.group(1)) if mr else -1,
        'validation': status,
        '_hlahd_cmd': hl_cmd,
        '_extract_cmd': mask(files.get('task/extract.command.sh', ''), sid),
    }


def main():
    arch = C.D_INPUTS / 'typing'
    man = pd.read_csv(C.D_INPUTS / 'archive_manifest.tsv', sep='\t', dtype=str)
    ids = [ln.split()[1] for ln in C.COHORT_KEEP.read_text().splitlines() if ln.strip()]
    if sorted(man.sample_id) != sorted(ids):
        raise SystemExit(f'ABORT: archive manifest has {len(man)} samples, cohort {len(ids)}; '
                         f'they differ')
    paths = [arch / f'{s}.tar.gz' for s in ids]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(one, paths, chunksize=16))
    df = pd.DataFrame(rows)

    vers = Counter(df.hlahd_version)
    hl_cmds = Counter(df._hlahd_cmd)
    ex_cmds = Counter(df._extract_cmd)
    bad = []
    if len(vers) != 1 or '' in vers:
        bad.append(f'HLA-HD versions: {dict(vers)}')
    if len(hl_cmds) != 1:
        bad.append(f'{len(hl_cmds)} distinct HLA-HD commands')
    if len(ex_cmds) != 1:
        bad.append(f'{len(ex_cmds)} distinct extraction commands')
    if (df.validation != 'ok').any():
        bad.append(f'validation not ok: {int((df.validation != "ok").sum())} samples')
    if (df.read_pairs < 0).any():
        bad.append(f'no read count for {int((df.read_pairs < 0).sum())} samples')
    if bad:
        raise SystemExit('ABORT: the earlier run is not one consistent typing: ' + '; '.join(bad))

    meta = pd.read_csv(C.OLD_MANIFEST, sep='\t', dtype=str)
    meta = meta[['sample_id', 'group', 'sex', 'platform', 'observed_depth']]
    out = df.drop(columns=['_hlahd_cmd', '_extract_cmd']).merge(meta, on='sample_id', how='left')
    if out.group.isna().any():
        raise SystemExit('ABORT: samples without metadata in the earlier manifest')
    out = out.merge(man[['sample_id', 'md5']].rename(columns={'md5': 'archive_md5'}), on='sample_id')
    out.to_csv(C.D_INPUTS / 'sample_inputs.tsv', sep='\t', index=False)

    def tomd5(sumfile, name):
        for ln in Path(sumfile).read_text().splitlines():
            f = ln.split()
            if len(f) == 2 and f[1].endswith(name):
                return f[0]
        return ''

    pinned = {
        'hlahd_freq_data': {'path': str(C.HLAHD_FREQ), 'md5_dir': md5_dir(C.HLAHD_FREQ)},
        'hlahd_gene_split': {'path': str(C.HLAHD_SPLIT), 'md5': md5_file(C.HLAHD_SPLIT)},
        'hlahd_dictionary': {'path': str(C.HLAHD_DICT)},
        'imgt_alignments': {'path': str(C.IMGT_ALIGN), 'md5_dir': md5_dir(C.IMGT_ALIGN)},
        'p_groups': {'path': str(C.P_GROUPS), 'md5': md5_file(C.P_GROUPS)},
        'imgt_allele_history': {'path': str(C.IMGT_HISTORY), 'md5': md5_file(C.IMGT_HISTORY),
                                'source': 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/3640/Allelelist_history.txt'},
        'tommo_panel': {'path': str(C.TOMMO), 'md5': md5_file(C.TOMMO),
                        'md5_as_downloaded': tomd5(C.TOMMO_MD5, C.TOMMO.name)},
        'kg_panel': {'path': str(C.KG), 'md5': md5_file(C.KG),
                     'md5_as_downloaded': tomd5(C.KG_MD5, C.KG.name)},
        'cohort_keep': {'path': str(C.COHORT_KEEP), 'md5': md5_file(C.COHORT_KEEP)},
        'earlier_manifest': {'path': str(C.OLD_MANIFEST), 'md5': md5_file(C.OLD_MANIFEST)},
        'region_list': {'path': str(C.OLD_CONTIGS), 'md5': md5_file(C.OLD_CONTIGS)},
    }
    for k in ('tommo_panel', 'kg_panel'):
        if pinned[k]['md5_as_downloaded'] and pinned[k]['md5'] != pinned[k]['md5_as_downloaded']:
            raise SystemExit(f'ABORT: {k} differs from its recorded download md5')
    # HLA-HD's tie-breaking prior is pinned by copy: it lives in the installation
    # directory, which this component does not own.
    pin = C.D_INPUTS / 'hlahd_freq_data'
    if not pin.exists():
        shutil.copytree(C.HLAHD_FREQ, pin)
    if md5_dir(pin) != pinned['hlahd_freq_data']['md5_dir']:
        raise SystemExit(f'ABORT: {pin} differs from {C.HLAHD_FREQ}')
    pinned['hlahd_freq_data']['pinned_copy'] = str(pin)
    prov = {
        'samples': len(out),
        'groups': out.group.value_counts().to_dict(),
        'hlahd_version': next(iter(vers)),
        'hlahd_command': next(iter(hl_cmds)),
        'extract_command': next(iter(ex_cmds)),
        'regions': int(out.regions.iloc[0]),
        'read_pairs': {'median': int(out.read_pairs.median()),
                       'p1': int(out.read_pairs.quantile(0.01)),
                       'min': int(out.read_pairs.min()), 'max': int(out.read_pairs.max())},
        'validation_ok': int((out.validation == 'ok').sum()),
        'pinned': pinned,
    }
    (C.D_INPUTS / 'provenance.json').write_text(json.dumps(prov, indent=2) + '\n')
    print(json.dumps({k: v for k, v in prov.items() if k != 'pinned'}, indent=2))


if __name__ == '__main__':
    main()
