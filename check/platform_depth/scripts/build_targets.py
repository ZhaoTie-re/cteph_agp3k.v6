#!/usr/bin/env python3
"""Stage 02 — what is re-called and re-tested.

Signal variants
  GWAS   the 168 signal variants of the 16 common-variant loci (leads, alternative leads,
         credible-set members; integration/saige_signals 00.selection/variants.tsv).
  genes  the Tier 1-2 genes of assoc_gene (robust_genes/tiers.tsv) in BOTH variant sets:
         per cohort x stratum, exactly the variants rvtest tested — the positions on the gene's
         set-file line that are in the cohort's tested-variant list.

Calling intervals
  Each signal variant +/- FLANK, merged, clipped to the original calling intervals and aligned
  to their 5-kb read-shard grid (HaplotypeCaller's results depend on interval boundaries; on the
  grid the shards are the original run's shards).

Fixed sites
  Every variant the analyses read inside the calling intervals: random-model variants of any
  cohort (the GWAS signals and their regional context) and the genes' tested variants. Re-called
  genotypes are kept at these IDs only; the variant set never changes.

Null SNPs (platform analysis)
  Random LD-pruned common SNVs and indels (assoc_saige's GRM markers, present in every cohort),
  outside the MHC, the 16 loci and the calling intervals, drawn per MAF bin: the reference
  distribution for Q and the check of the shared-control correlation.
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402

STRATA = ['moderate_high', 'low_moderate_high']
MHC = ('chr6', 28_510_120, 33_480_577)          # GRCh38 extended MHC, as assoc_hla.v2
MAF_BINS = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50]


def gwas_signals():
    v = K.read_tsv(K.SAIGE_SIGNALS / '00.selection' / 'variants.tsv')
    if len(v) != 168 or v.variant_id.duplicated().any():
        sys.exit(f'ABORT: expected 168 distinct GWAS signal variants, found {len(v)}')
    loci = K.read_tsv(K.SAIGE_SIGNALS / '00.selection' / 'loci.tsv')
    v = v.merge(loci[['locus_id', 'label', 'tier']].rename(columns={'label': 'locus', 'tier': 'locus_tier'}),
                on='locus_id', how='left')
    return pd.DataFrame({'kind': 'gwas', 'unit': v.locus, 'unit_id': v.locus_id, 'role': v.role,
                         'variant_id': v.variant_id, 'cohort': 'all', 'stratum': '', 'tier': v.locus_tier})


def gene_signals():
    tiers = K.read_tsv(K.ASSOC_GENE / '_comparison' / 'robust_genes' / 'tiers.tsv')
    t12 = tiers[tiers.tier.isin(['1', '2'])]
    genes = sorted(set(t12.gene_symbol))
    tier_of = {(r.gene_symbol, r.stratum): r.tier for r in tiers.itertuples()}
    chrom_of = dict(zip(t12.gene_symbol, t12.chrom))
    rows = []
    for c in K.COHORTS:
        for s in STRATA:
            tested = set((K.ASSOC_GENE / c / '01.map' / f'variants.{s}.txt').read_text().split())
            by_pos = defaultdict(list)
            for vid in tested:
                ch, p, _, _ = K.split_id(vid)
                by_pos[(ch, p)].append(vid)
            for g in genes:
                setf = K.ASSOC_GENE / c / '02.test_inputs' / s / f'rvtest.set.chr{chrom_of[g]}.txt'
                line = [ln for ln in setf.read_text().splitlines() if ln.split('\t')[0] == g]
                if len(line) != 1:
                    sys.exit(f'ABORT: {g} has {len(line)} lines in {setf}')
                for rng in line[0].split('\t')[1].split(','):
                    ch, span = rng.split(':')
                    a, b = (int(x) for x in span.split('-'))
                    if a != b:
                        sys.exit(f'ABORT: {g} set range {rng} is not a single position')
                    ids = by_pos.get((f'chr{ch}', a), [])
                    if not ids:
                        sys.exit(f'ABORT: {g} {c} {s}: no tested variant at chr{ch}:{a}')
                    for vid in ids:
                        rows.append({'kind': 'gene', 'unit': g, 'unit_id': g, 'role': 'tested',
                                     'variant_id': vid, 'cohort': c, 'stratum': s,
                                     'tier': tier_of.get((g, s), '0')})
    return pd.DataFrame(rows), genes


def calling_intervals(chroms):
    iv = {}
    for ch in chroms:
        rows = []
        for ln in (K.CALLING_INTERVALS / f'wgs_calling_regions.hg38.{ch}.interval_list').read_text().splitlines():
            if ln.startswith('@') or not ln.strip():
                continue
            f = ln.split('\t')
            rows.append((int(f[1]), int(f[2])))
        iv[ch] = sorted(rows)
    return iv


def merge(spans):
    out = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def hc_intervals(positions, flank):
    """positions: chrom -> [pos]; returns chrom -> [(start, end, calling_start)], 1-based inclusive."""
    civ = calling_intervals(sorted(positions))
    out = {}
    for ch, ps in positions.items():
        raw = merge([(p - flank, p + flank) for p in ps])
        pieces = []
        for s, e in raw:
            for s0, e0 in civ[ch]:
                lo, hi = max(s, s0), min(e, e0)
                if lo > hi:
                    continue
                a = s0 + ((lo - s0) // K.SHARD) * K.SHARD
                b = min(e0, s0 + (-(-(hi - s0 + 1) // K.SHARD)) * K.SHARD - 1)
                pieces.append((a, b, s0))
        # merge within the same calling interval only, so every piece stays on its own grid
        merged = []
        for a, b, s0 in sorted(pieces):
            if merged and merged[-1][2] == s0 and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b), s0)
            else:
                merged.append((a, b, s0))
        for p in ps:
            if not any(a <= p <= b for a, b, _ in merged):
                sys.exit(f'ABORT: {ch}:{p} falls outside the calling intervals')
        out[ch] = merged
    return out


def chrom_key(ch):
    return int(ch.replace('chr', ''))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--flank', type=int, default=10_000)
    ap.add_argument('--extract-pad', type=int, default=1_000)
    ap.add_argument('--null-per-bin', type=int, default=3000, help='SNVs per MAF bin')
    ap.add_argument('--null-indel-per-bin', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=20261030)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    gw = gwas_signals()
    gn, genes = gene_signals()
    sig = pd.concat([gw, gn], ignore_index=True)
    sig[['chrom', 'pos', 'ref', 'alt']] = pd.DataFrame([K.split_id(v) for v in sig.variant_id], index=sig.index)
    K.write_tsv(sig, out / 'signal_variants.tsv')

    positions = defaultdict(set)
    for ch, p in zip(sig.chrom, sig.pos):
        positions[ch].add(int(p))
    hc = hc_intervals({k: sorted(v) for k, v in positions.items()}, a.flank)

    with open(out / 'hc.intervals', 'w') as fh:          # GATK -L, 1-based inclusive
        for ch in sorted(hc, key=chrom_key):
            for s, e, _ in hc[ch]:
                fh.write(f'{ch}:{s}-{e}\n')
    iv = pd.DataFrame([{'chrom': ch, 'start': s, 'end': e, 'calling_interval_start': s0,
                        'grid_offset_start': (s - s0) % K.SHARD, 'length': e - s + 1}
                       for ch in sorted(hc, key=chrom_key) for s, e, s0 in hc[ch]])
    K.write_tsv(iv, out / 'hc_intervals.tsv')
    with open(out / 'extract.bed', 'w') as fh:            # samtools --region-file, 0-based
        for ch in sorted(hc, key=chrom_key):
            for s, e in merge([(s - a.extract_pad, e + a.extract_pad) for s, e, _ in hc[ch]]):
                fh.write(f'{ch}\t{max(0, s - 1)}\t{e}\n')

    # fixed sites: random-model variants of any cohort inside the intervals + the genes' tested variants
    def inside(ch, p):
        return any(s <= p <= e for s, e, _ in hc.get(ch, []))

    fixed = {}
    member = defaultdict(dict)
    for c in K.COHORTS:
        bim = pd.read_csv(K.MODEL_INPUTS / c / 'genotype' / 'random_model' / f'{K.RANDOM_PREFIX}.bim',
                          sep='\t', header=None, dtype=str, usecols=[0, 1, 3])
        bim.columns = ['chr', 'id', 'pos']
        bim = bim[bim.chr.isin([ch.replace('chr', '') for ch in hc])]
        for vid, ch, p in zip(bim.id, bim.chr, bim.pos.astype(int)):
            if inside('chr' + ch, p):
                fixed.setdefault(vid, 'gwas_context')
                member[vid][f'random_{c}'] = 1
        rvb = pd.read_csv(K.TUNING_RV / c / '02.callset_filter' / 'filtered.bim', sep='\t', header=None,
                          dtype=str, usecols=[1])[1]
        gene_ids = set(gn.variant_id[gn.cohort == c])
        missing = gene_ids - set(rvb)
        if missing:
            sys.exit(f'ABORT: {len(missing)} tested gene variants are not in {c} filtered.bim')
        for vid in gene_ids:
            member[vid][f'rv_{c}'] = 1
    for vid in gw.variant_id:
        if vid not in fixed:
            sys.exit(f'ABORT: GWAS signal {vid} is in no random-model fileset')
        fixed[vid] = 'gwas_signal'
    for vid in gn.variant_id.unique():
        fixed[vid] = 'gene' if vid not in fixed else fixed[vid] + ';gene'
    fs = []
    for vid, kind in fixed.items():
        ch, p, r, al = K.split_id(vid)
        row = {'variant_id': vid, 'chrom': ch, 'pos': p, 'ref': r, 'alt': al, 'kind': kind}
        for c in K.COHORTS:
            row[f'random_{c}'] = member[vid].get(f'random_{c}', 0)
            row[f'rv_{c}'] = member[vid].get(f'rv_{c}', 0)
        fs.append(row)
    fs = pd.DataFrame(fs)
    fs['_k'] = fs.chrom.map(chrom_key)
    fs = fs.sort_values(['_k', 'pos', 'variant_id']).drop(columns='_k')
    K.write_tsv(fs, out / 'fixed_sites.tsv')
    with open(out / 'fixed_positions.intervals', 'w') as fh:   # GenotypeGVCFs -L, one base each
        for ch, p in fs[['chrom', 'pos']].drop_duplicates().itertuples(index=False):
            fh.write(f'{ch}:{p}-{p}\n')

    # null SNPs for the platform analysis
    prune = None
    for c in K.COHORTS:
        s = set((K.ASSOC_SAIGE / c / '00.prep' / 'grm.prune.in').read_text().split())
        prune = s if prune is None else prune & s
    af = pd.read_csv(K.ASSOC_SAIGE / 'full_mainland' / '01.assoc' / 'main' / 'saige.assoc.tsv', sep='\t',
                     usecols=['CHR', 'POS', 'MarkerID', 'AF_Allele2'], dtype={'MarkerID': str})
    af = af[af.MarkerID.isin(prune)].copy()
    af['chrom'] = 'chr' + af.CHR.astype(str)
    loci = K.read_tsv(K.SAIGE_SIGNALS / '00.selection' / 'loci.tsv')
    excl = [MHC]
    for vid in loci.variant_id:
        ch, p, _, _ = K.split_id(vid)
        excl.append((ch, p - 1_000_000, p + 1_000_000))
    for ch in hc:
        excl += [(ch, s, e) for s, e, _ in hc[ch]]
    keep = np.ones(len(af), bool)
    for ch, s, e in excl:
        keep &= ~((af.chrom.values == ch) & (af.POS.values >= s) & (af.POS.values <= e))
    af = af[keep].copy()
    af['maf'] = np.minimum(af.AF_Allele2, 1 - af.AF_Allele2)
    ra = af.MarkerID.str.split(':', expand=True)
    af['is_indel'] = (ra[2].str.len() != 1) | (ra[3].str.len() != 1)
    af['maf_bin'] = pd.cut(af.maf, MAF_BINS, right=False, labels=[f'{lo}-{hi}' for lo, hi in
                                                                    zip(MAF_BINS[:-1], MAF_BINS[1:])])
    rng = np.random.default_rng(a.seed)
    picks = []
    for (b, ind), g in af.dropna(subset=['maf_bin']).groupby(['maf_bin', 'is_indel'], observed=True):
        n = a.null_indel_per_bin if ind else a.null_per_bin
        idx = rng.choice(len(g), size=min(n, len(g)), replace=False)
        picks.append(g.iloc[np.sort(idx)])
    nul = pd.concat(picks)[['MarkerID', 'chrom', 'POS', 'maf', 'maf_bin', 'is_indel']]
    nul.columns = ['variant_id', 'chrom', 'pos', 'maf_full', 'maf_bin', 'is_indel']
    K.write_tsv(nul.sort_values(['chrom', 'pos']), out / 'null_snps.tsv')

    qc = {
        'gwas_signal_variants': int((sig.kind == 'gwas').sum()),
        'genes': genes,
        'gene_tested_variants_distinct': int(gn.variant_id.nunique()),
        'gene_rows': gn.groupby(['unit', 'stratum']).variant_id.nunique().reset_index()
                       .apply(lambda r: f'{r.unit}/{r.stratum}={r.variant_id}', axis=1).tolist(),
        'hc_intervals': len(iv), 'hc_bp': int(iv.length.sum()),
        'fixed_sites': len(fs), 'fixed_by_kind': fs.kind.value_counts().to_dict(),
        'null_snps': len(nul),
        'null_by_bin': [f'{b}/{"indel" if i else "snv"}={n}' for (b, i), n in
                        nul.groupby(['maf_bin', 'is_indel'], observed=True).size().items()],
        'flank': a.flank, 'shard': K.SHARD, 'seed': a.seed,
    }
    (out / 'targets_qc.json').write_text(json.dumps(qc, indent=2) + '\n')
    print(json.dumps(qc, indent=2))


if __name__ == '__main__':
    main()
