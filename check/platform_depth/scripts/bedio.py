#!/usr/bin/env python3
"""PLINK 1 .bed rows, read and written byte for byte (no plink round trip, so allele order and
every untouched genotype stay exactly as the analysis read them).

  subset   copy the rows of the given variant IDs (in .bim order) into a new fileset
  swap     replace some samples' genotypes with those of a re-called VCF, matched by variant ID
  compare  genotype agreement of two filesets with the same variants, for some samples

.bed codes per sample (2 bits, SNP-major): 00 = homozygous for the .bim column-5 allele,
01 = missing, 10 = heterozygous, 11 = homozygous for the column-6 allele.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

MAGIC = bytes([0x6C, 0x1B, 0x01])
HOM5, MISS, HET, HOM6 = 0, 1, 2, 3


def read_fam(prefix):
    return [ln.split()[1] for ln in Path(f'{prefix}.fam').read_text().splitlines() if ln.strip()]


def read_bim(prefix):
    rows = []
    with open(f'{prefix}.bim') as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            rows.append(f)
    return rows


def bpv(n):
    return (n + 3) // 4


def read_rows(prefix, idx, n):
    """Raw bytes of the given variant rows (sorted indices)."""
    b = bpv(n)
    out = np.empty((len(idx), b), np.uint8)
    with open(f'{prefix}.bed', 'rb') as fh:
        if fh.read(3) != MAGIC:
            sys.exit(f'ABORT: {prefix}.bed is not a SNP-major PLINK 1 bed')
        for k, i in enumerate(idx):
            fh.seek(3 + i * b)
            out[k] = np.frombuffer(fh.read(b), np.uint8)
    return out


def unpack(rows, n):
    """bytes (m, bpv) -> codes (m, n)."""
    m = rows.shape[0]
    codes = np.empty((m, rows.shape[1] * 4), np.uint8)
    for s in range(4):
        codes[:, s::4] = (rows >> (2 * s)) & 3
    return codes[:, :n]


def pack(codes):
    m, n = codes.shape
    b = bpv(n)
    pad = np.zeros((m, b * 4), np.uint8)
    pad[:, :n] = codes
    out = np.zeros((m, b), np.uint8)
    for s in range(4):
        out |= (pad[:, s::4] << (2 * s)).astype(np.uint8)
    return out


def write_fileset(prefix, rows, bim_rows, fam_src):
    with open(f'{prefix}.bed', 'wb') as fh:
        fh.write(MAGIC)
        fh.write(rows.tobytes())
    with open(f'{prefix}.bim', 'w') as fh:
        for f in bim_rows:
            fh.write('\t'.join(f) + '\n')
    Path(f'{prefix}.fam').write_text(Path(f'{fam_src}.fam').read_text())


def cmd_subset(a):
    want = set(Path(a.ids).read_text().split())
    bim = read_bim(a.bfile)
    idx = [i for i, f in enumerate(bim) if f[1] in want]
    found = {bim[i][1] for i in idx}
    if a.require_all and found != want:
        sys.exit(f'ABORT: {len(want - found)} requested IDs are not in {a.bfile}.bim')
    n = len(read_fam(a.bfile))
    rows = read_rows(a.bfile, idx, n)
    write_fileset(a.out, rows, [bim[i] for i in idx], a.bfile)
    print(f'[subset] {len(idx)} of {len(want)} IDs from {a.bfile} ({n} samples) -> {a.out}')


def alt_counts(vcf, samples, ids):
    """{variant_id: np.array of ALT counts (-1 = missing)} for the given VCF samples."""
    import pysam
    vf = pysam.VariantFile(vcf)
    names = list(vf.header.samples)
    miss = [s for s in samples if s not in names]
    if miss:
        sys.exit(f'ABORT: {len(miss)} samples are not in {vcf}, e.g. {miss[:3]}')
    out = {}
    for rec in vf.fetch():
        if rec.id not in ids:
            continue
        if len(rec.alts or ()) != 1:
            sys.exit(f'ABORT: {rec.id} is not biallelic in {vcf}')
        if rec.id in out:
            sys.exit(f'ABORT: {rec.id} appears twice in {vcf}')
        g = np.full(len(samples), -1, np.int8)
        for k, s in enumerate(samples):
            gt = rec.samples[s]['GT']
            if gt is None or any(x is None for x in gt):
                continue
            g[k] = sum(1 for x in gt if x == 1)
            if any(x not in (0, 1) for x in gt):
                sys.exit(f'ABORT: {rec.id} {s} has genotype {gt}')
        out[rec.id] = g
    return out


def cmd_swap(a):
    fam = read_fam(a.bfile)
    n = len(fam)
    bim = read_bim(a.bfile)
    pairs = [ln.split('\t') for ln in Path(a.samples).read_text().splitlines() if ln.strip()]
    sample_col = {vcf_name: fam.index(sid) for sid, vcf_name in pairs if sid in fam}
    if not sample_col:
        sys.exit('ABORT: none of the samples to swap is in the fileset')
    vcf_names = list(sample_col)
    ids = {f[1] for f in bim}
    calls = alt_counts(a.vcf, vcf_names, ids)
    missing = ids - set(calls)
    if missing:
        sys.exit(f'ABORT: {len(missing)} variants have no record in {a.vcf}, e.g. {sorted(missing)[:3]}')
    rows = read_rows(a.bfile, list(range(len(bim))), n)
    codes = unpack(rows, n)
    before = codes.copy()
    cols = np.array([sample_col[v] for v in vcf_names])
    for i, f in enumerate(bim):
        _, _, ref, alt = f[1].split(':', 3)
        a5, a6 = f[4], f[5]
        g = calls[f[1]]
        if (a5, a6) == (alt, ref):
            lut = {2: HOM5, 1: HET, 0: HOM6, -1: MISS}
        elif (a5, a6) == (ref, alt):
            lut = {0: HOM5, 1: HET, 2: HOM6, -1: MISS}
        else:
            sys.exit(f'ABORT: {f[1]} has .bim alleles {a5}/{a6}')
        codes[i, cols] = np.vectorize(lut.get, otypes=[np.uint8])(g)
    untouched = np.setdiff1d(np.arange(n), cols)
    if not np.array_equal(codes[:, untouched], before[:, untouched]):
        sys.exit('ABORT: the swap changed an untouched sample')
    write_fileset(a.out, pack(codes), bim, a.bfile)
    changed = (codes[:, cols] != before[:, cols])
    rep = {'variants': len(bim), 'samples_swapped': int(len(cols)), 'samples_untouched': int(len(untouched)),
           'genotypes_changed': int(changed.sum()), 'share_changed': float(changed.mean()),
           'untouched_identical': True}
    if a.report:
        Path(a.report).write_text(json.dumps(rep, indent=2) + '\n')
    print(f'[swap] {json.dumps(rep)}')


def cmd_compare(a):
    fa, fb = read_fam(a.a), read_fam(a.b)
    ba, bb = read_bim(a.a), read_bim(a.b)
    if [f[1] for f in ba] != [f[1] for f in bb] or [f[4:6] for f in ba] != [f[4:6] for f in bb]:
        sys.exit('ABORT: the two filesets differ in variants or allele order')
    want = [s for s in Path(a.samples).read_text().split() if s in fa and s in fb]
    ca = unpack(read_rows(a.a, list(range(len(ba))), len(fa)), len(fa))[:, [fa.index(s) for s in want]]
    cb = unpack(read_rows(a.b, list(range(len(bb))), len(fb)), len(fb))[:, [fb.index(s) for s in want]]
    labels = {HOM5: 'hom5', MISS: 'miss', HET: 'het', HOM6: 'hom6'}
    with open(a.out, 'w') as fh:
        fh.write('variant_id\tn\tsame\tdiff_called\tmiss_a\tmiss_b\tlost\tgained\n')
        for i, f in enumerate(ba):
            x, y = ca[i], cb[i]
            same = int((x == y).sum())
            lost = int(((x != MISS) & (y == MISS)).sum())
            gained = int(((x == MISS) & (y != MISS)).sum())
            diff = int(((x != MISS) & (y != MISS) & (x != y)).sum())
            fh.write(f'{f[1]}\t{len(want)}\t{same}\t{diff}\t{int((x == MISS).sum())}\t{int((y == MISS).sum())}\t'
                     f'{lost}\t{gained}\n')
    print(f'[compare] {len(ba)} variants x {len(want)} samples -> {a.out} ({labels})')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('subset')
    p.add_argument('--bfile', required=True)
    p.add_argument('--ids', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--require-all', action='store_true')
    p = sub.add_parser('swap')
    p.add_argument('--bfile', required=True)
    p.add_argument('--vcf', required=True)
    p.add_argument('--samples', required=True, help='TSV: sample_id <TAB> name in the VCF')
    p.add_argument('--out', required=True)
    p.add_argument('--report')
    p = sub.add_parser('compare')
    p.add_argument('--a', required=True)
    p.add_argument('--b', required=True)
    p.add_argument('--samples', required=True)
    p.add_argument('--out', required=True)
    a = ap.parse_args()
    {'subset': cmd_subset, 'swap': cmd_swap, 'compare': cmd_compare}[a.cmd](a)


if __name__ == '__main__':
    main()
