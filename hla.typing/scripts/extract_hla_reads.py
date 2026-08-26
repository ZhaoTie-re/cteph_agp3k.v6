#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One CRAM -> the pair of FASTQ files HLA-HD reads.
#
#           Three things here are not optional:
#
#           1. The index is given explicitly and passed with `samtools view -X`.
#              Deriving it as <cram>.crai fails for 3,117 of this study's samples,
#              whose CRAM path is a symlink and whose index lives beside the real
#              file. Nothing is ever written next to the source CRAM.
#
#           2. `samtools collate` before `samtools fastq`. The extracted BAM is
#              coordinate-sorted, so without collating, R1 and R2 come out in
#              different orders AND different counts. HLA-HD survives that because
#              it aligns the two files separately and re-pairs by read name, but
#              the files are then wrong for every other consumer, and silently so.
#
#           3. Every exit code is checked, INCLUDING the decompressor's. The read
#              counts that prove R1 and R2 are paired were once taken from two bare
#              `Popen(['zcat', …])` whose return code was never read: if zcat failed,
#              both counts came out 0, `n1 != n2` was False, and an empty FASTQ pair
#              passed the check that exists to catch exactly that.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--cram', required=True)
    p.add_argument('--crai', required=True, help='resolved index; passed to samtools -X')
    p.add_argument('--reference', required=True)
    p.add_argument('--contig-list', required=True)
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--out-prefix', default=None, help='default: <sample>.hla')
    return p.parse_args()


def run(cmd, **kw):
    """Run, and fail loudly with the command that failed."""
    r = subprocess.run([str(c) for c in cmd], **kw)
    if r.returncode != 0:
        raise SystemExit(f'ABORT: exit {r.returncode} from: {" ".join(str(c) for c in cmd)}')
    return r


def count_reads(path):
    """Read pairs in a gzipped FASTQ, with the decompressor's exit code checked.

    An unchecked decompressor returns 0 lines on failure, which is indistinguishable
    from an empty file and passes every equality test downstream.
    """
    z = subprocess.Popen(['gzip', '-cd', str(path)], stdout=subprocess.PIPE)
    n = sum(1 for _ in z.stdout)
    z.stdout.close()
    rc = z.wait()
    if rc != 0:
        raise SystemExit(f'ABORT: gzip -cd exited {rc} on {path}')
    if n % 4:
        raise SystemExit(f'ABORT: {path} has {n} lines, not a multiple of 4 — truncated')
    return n // 4


def main():
    args = parse_args()
    prefix = args.out_prefix or f'{args.sample}.hla'
    for p in (args.cram, args.crai, args.reference, args.contig_list):
        if not Path(p).exists():
            raise SystemExit(f'ABORT: missing input {p}')

    regions = [ln.strip() for ln in Path(args.contig_list).read_text().splitlines() if ln.strip()]
    if not regions:
        raise SystemExit(f'ABORT: {args.contig_list} is empty')

    bam = f'{prefix}.regions.bam'
    r1, r2 = f'{prefix}.R1.fastq.gz', f'{prefix}.R2.fastq.gz'

    # -X makes the explicit index argument legal; it is the whole reason this works
    # on a symlinked CRAM whose index is elsewhere.
    run([args.samtools, 'view', '-T', args.reference, '-@', args.threads,
         '-b', '-X', args.cram, args.crai, '-o', bam] + regions)

    # -O streams collated groups to stdout without a temp-file shuffle.
    collate = subprocess.Popen(
        [args.samtools, 'collate', '-u', '-O', '-@', str(args.threads), bam, f'{prefix}.collate.tmp'],
        stdout=subprocess.PIPE)
    fastq = subprocess.Popen(
        [args.samtools, 'fastq', '-@', str(args.threads),
         '-1', r1, '-2', r2,
         '-0', '/dev/null',          # neither/both READ1+READ2 set
         '-s', f'{prefix}.singleton.fastq.gz',   # named, not leaked to stdout
         '-n',                       # keep the name as-is; no /1 /2 suffix to undo
         '-'],
        stdin=collate.stdout)
    collate.stdout.close()
    fq_rc = fastq.wait()
    col_rc = collate.wait()
    if col_rc != 0:
        raise SystemExit(f'ABORT: samtools collate exited {col_rc}')
    if fq_rc != 0:
        raise SystemExit(f'ABORT: samtools fastq exited {fq_rc}')

    for f in (r1, r2):
        if not Path(f).exists() or Path(f).stat().st_size < 100:
            raise SystemExit(f'ABORT: {f} is missing or empty — extraction produced no reads')

    Path(bam).unlink(missing_ok=True)
    for tmp in Path('.').glob(f'{prefix}.collate.tmp*'):
        tmp.unlink(missing_ok=True)

    n1, n2 = count_reads(r1), count_reads(r2)
    if n1 != n2:
        raise SystemExit(f'ABORT: {args.sample} R1 has {n1} reads and R2 has {n2}. '
                         f'Collating should have made these equal.')
    print(f'[extract_hla_reads] {args.sample}: {n1:,} read pairs over '
          f'{len(regions)} region(s) -> {r1}, {r2}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in extract_hla_reads: {e}', file=sys.stderr)
        sys.exit(1)
