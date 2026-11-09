#!/usr/bin/env python3
"""Step 00a — archive the earlier HLA-HD output, one tar.gz per sample.

The typing itself is NOT re-run. This step copies, for every cohort sample, the
files the later steps read and the task records that show how they were made:

    result/<s>_final.result.txt       HLA-HD's reported genotype per gene
    result/<s>_<gene>.est.txt         tied best allele pairs and their exon scores
    result/<s>_<gene>.read.txt        reads assigned to each allele
    hlahd.run.log, typing_validation.tsv
    task/hlahd.command.sh|.trace      the exact HLA-HD command and its run record
    task/extract.command.sh|.out|.err|.trace
                                      the read-extraction command and its read counts

The source directories are symlinks into hla.typing/work/, which also holds
~0.8 GB of intermediate alignments per sample that are not needed here. One
archive per sample (about 0.5 MB) replaces ~70 small files on Lustre and makes the
inputs of this component independent of that work/ tree.

Writes results/00.inputs/typing/<s>.tar.gz and results/00.inputs/archive_manifest.tsv
(sample, md5, members, source paths). Idempotent: an archive that already exists is
not rebuilt (archives are written to a temporary name and renamed, so an existing
one is complete); the manifest is always rebuilt from the archives on disk.
"""
import argparse
import hashlib
import io
import os
import sys
import tarfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

TASK_FILES = ('.command.sh', '.command.trace')
EXTRACT_FILES = ('.command.sh', '.command.out', '.command.err', '.command.trace')


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def cohort_ids():
    ids = []
    for ln in C.COHORT_KEEP.read_text().splitlines():
        f = ln.split()
        if f:
            ids.append(f[1] if len(f) > 1 else f[0])
    if len(ids) != len(set(ids)):
        raise SystemExit(f'ABORT: duplicate IDs in {C.COHORT_KEEP}')
    return ids


def archive_one(sid, out_dir):
    """Build <sid>.tar.gz; return (sid, md5, n_members, typing_src, extract_src)."""
    typ = C.OLD_TYPING / sid
    if not typ.exists():
        raise RuntimeError(f'{sid}: no typing directory {typ}')
    typ_real = typ.resolve()
    hlahd_task = typ_real.parent
    r1 = C.OLD_READS / f'{sid}.hla.R1.fastq.gz'
    if not r1.exists():
        raise RuntimeError(f'{sid}: no extracted reads {r1}')
    ext_task = r1.resolve().parent

    res = typ_real / 'result'
    members = []
    final = res / f'{sid}_final.result.txt'
    if not final.exists():
        raise RuntimeError(f'{sid}: missing {final}')
    members.append((final, f'result/{final.name}'))
    for p in sorted(res.glob(f'{sid}_*.est.txt')) + sorted(res.glob(f'{sid}_*.read.txt')):
        members.append((p, f'result/{p.name}'))
    for name in ('hlahd.run.log', 'typing_validation.tsv'):
        p = typ_real / name
        if not p.exists():
            raise RuntimeError(f'{sid}: missing {p}')
        members.append((p, name))
    for name in TASK_FILES:
        p = hlahd_task / name
        if p.exists():
            members.append((p, f'task/hlahd{name}'))
    for name in EXTRACT_FILES:
        p = ext_task / name
        if p.exists():
            members.append((p, f'task/extract{name}'))

    out = out_dir / f'{sid}.tar.gz'
    tmp = out_dir / f'.{sid}.tar.gz.tmp'
    # Fixed mtime/owner so the archive's bytes depend on content only.
    with tarfile.open(tmp, 'w:gz', compresslevel=6) as tf:
        for src, arc in members:
            data = src.read_bytes()
            ti = tarfile.TarInfo(name=f'{sid}/{arc}')
            ti.size = len(data)
            ti.mtime = int(src.stat().st_mtime)
            ti.mode = 0o644
            tf.addfile(ti, io.BytesIO(data))
    os.replace(tmp, out)
    return sid, md5(out), len(members), str(typ_real), str(ext_task)


def describe(sid, out_dir):
    """Manifest row for an archive already on disk."""
    p = out_dir / f'{sid}.tar.gz'
    with tarfile.open(p, 'r:gz') as tf:
        n = sum(1 for m in tf.getmembers() if m.isfile())
    typ = (C.OLD_TYPING / sid).resolve()
    ext = (C.OLD_READS / f'{sid}.hla.R1.fastq.gz').resolve().parent
    return sid, md5(p), n, str(typ), str(ext)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--limit', type=int, default=0, help='first N samples only (testing)')
    args = ap.parse_args()

    out_dir = C.D_INPUTS / 'typing'
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = C.D_INPUTS / 'archive_manifest.tsv'

    ids = cohort_ids()
    if args.limit:
        ids = ids[:args.limit]
    for tmp in out_dir.glob('.*.tar.gz.tmp'):
        tmp.unlink()                      # left by an interrupted run
    todo = [s for s in ids if not (out_dir / f'{s}.tar.gz').exists()]
    print(f'[archive] {len(ids)} cohort samples; {len(ids) - len(todo)} already archived; '
          f'{len(todo)} to build with {args.workers} workers', flush=True)

    rows, errors = {}, []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(archive_one, s, out_dir): s for s in todo}
        for i, fu in enumerate(as_completed(futs), 1):
            s = futs[fu]
            try:
                sid, h, n, ts, es = fu.result()
                rows[sid] = [sid, h, str(n), ts, es]
            except Exception as e:  # collected, reported, and fatal below
                errors.append(f'{s}: {e}')
            if i % 250 == 0:
                print(f'[archive] {i}/{len(todo)}', flush=True)
        # archives built by an earlier run: describe them from disk
        futs = {ex.submit(describe, s, out_dir): s for s in ids if s not in rows
                and (out_dir / f'{s}.tar.gz').exists()}
        for fu in as_completed(futs):
            sid, h, n, ts, es = fu.result()
            rows[sid] = [sid, h, str(n), ts, es]

    with manifest.open('w') as fh:
        fh.write('sample_id\tmd5\tn_members\ttyping_source\textract_source\n')
        for s in ids:
            if s in rows:
                fh.write('\t'.join(rows[s]) + '\n')
    if errors:
        print(f'[archive] {len(errors)} sample(s) FAILED:', file=sys.stderr)
        for e in errors[:20]:
            print('   ', e, file=sys.stderr)
        raise SystemExit(1)
    n_members = sorted({r[2] for r in rows.values()})
    print(f'[archive] {len(rows)} archives -> {out_dir}; members per archive: {n_members}')


if __name__ == '__main__':
    main()
