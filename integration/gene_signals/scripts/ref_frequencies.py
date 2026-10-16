#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : ALT-allele frequencies of every table variant in external
#           references, one subcommand per source (adapted from
#           integration/saige_signals; gnomAD comes with the other ANNOVAR
#           annotations in annovar_annotations.py):
#             tommo   ToMMo 61KJPN genotype-frequency VCF (INFO/AF) and its ID
#                     field, cross-checked against the 60KJPN normalised VCF
#             kg      1000 Genomes 30x (NYGC) genotypes subset to the EAS panel
#                     and its JPT members (AF = AC/AN), cross-checked against
#                     plink2 --freq on the project EAS panel and INFO/AF_EAS_unrel
#                     (rare variants are often absent: the checks run on those present)
#           Reference records are split and left-normalised (bcftools norm)
#           before matching on chr:pos:REF:ALT. The references' tabix indexes are
#           older than their data; htslib warns on stderr and exits 0, so only
#           exit codes are checked.
# Inputs  : variants.tsv; the reference paths (templates use @@CHROM@@)
# Outputs : tommo.tsv | kg.tsv  one row per variant
# ---------------------------------------------------------------------------
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gs_common import chrom_plain, load_variants, log, num, read_tsv, run, write_tsv  # noqa: E402

BCFTOOLS = TABIX = FASTA = PLINK2 = None      # set from the command line in main()


def tpl(template, chrom):
    return template.replace('@@CHROM@@', str(chrom))


def bcf_norm_query(vcf, region, fmt, samples=None, keep_info='INFO/AC,INFO/AN,INFO/AF'):
    """view -r [-S] | annotate -x (all but keep_info, all FORMAT) | norm -m -any -f FASTA | query.

    Normalising every slice through one code path makes an un-trimmed or
    multi-allelic reference record comparable with our chr:pos:REF:ALT key.
    Every INFO field the query does not read is dropped BEFORE norm: the ToMMo
    genotype-frequency records carry ~37 kB of per-transcript and per-platform
    INFO, and splitting them made `norm` peak at 4.5 GB on one 8-kb slice
    (25 MB without them; identical output). With -S, view recomputes AC/AN for
    the subset before the FORMAT fields are dropped.
    """
    view = [BCFTOOLS, 'view', '-r', region, '-Ou', vcf]
    if samples:
        view[2:2] = ['-S', str(samples)]
    steps = [view,
             [BCFTOOLS, 'annotate', '-x', f'^{keep_info},FORMAT', '-Ou'],
             [BCFTOOLS, 'norm', '-m', '-any', '-f', FASTA, '-c', 'w', '-Ou'],
             [BCFTOOLS, 'query', '-f', fmt]]
    procs, prev = [], None
    for i, cmd in enumerate(steps):
        last = i == len(steps) - 1
        p = subprocess.Popen(cmd, stdin=prev, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE if last else subprocess.DEVNULL, text=last)
        if prev is not None:
            prev.close()
        prev = p.stdout
        procs.append(p)
    out, err = procs[-1].communicate()
    for p in procs[:-1]:
        p.wait()
    codes = [p.returncode for p in procs]
    if any(codes):
        raise RuntimeError(f'bcftools pipeline failed for {vcf} {region} (exit codes view/annotate/norm/query '
                           f'{codes}; a killed stage shows as a truncated stream): {err[-2000:]}')
    return [l.split('\t') for l in out.splitlines() if l]


def raw_alleles(vcf, chrom, pos):
    """(REF, ALT_k) pairs present in the UN-normalised records at a position."""
    pairs = set()
    for line in run([TABIX, vcf, f'chr{chrom}:{pos}-{pos}']).splitlines():
        f = line.split('\t')
        if len(f) >= 5 and f[1] == str(pos):
            for a in f[4].split(','):
                pairs.add((f[3], a))
    return pairs


def locus_windows(v):
    """One bcftools region per gene, spanning its variants +-100 bp."""
    w = {}
    for locus, g in v.groupby('locus_id'):
        chrom = g.chrom.iloc[0]
        assert g.chrom.nunique() == 1
        w[locus] = (chrom, f'chr{chrom}:{max(1, g.pos.min() - 100)}-{g.pos.max() + 100}')
    return w


def unique_variants(path):
    """variants.tsv -> one row per variant, its gene as the window key (a variant of two
    table genes is looked up once)."""
    v = load_variants(path).drop_duplicates('variant_id').copy()
    v['locus_id'] = v.set_name
    v['order'] = v.gene_order
    return v.reset_index(drop=True)


def fetch_tommo(v, template, crosscheck):
    rows = {}
    for locus, (chrom, region) in locus_windows(v).items():
        vcf = tpl(template, chrom)
        recs = bcf_norm_query(vcf, region,
                              '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/AN\t%INFO/AC\t%INFO/AF\n')
        by_key = {}
        for c, p, rsid, ref, alt, filt, an, ac, af in recs:
            by_key.setdefault(f'chr{chrom_plain(c)}:{p}:{ref}:{alt}', []).append((rsid, filt, an, ac, af))
        for _, r in v[v.locus_id == locus].iterrows():
            hits = by_key.get(r.variant_id, [])
            if not hits:
                rows[r.variant_id] = {'tommo_af': np.nan, 'tommo_an': np.nan, 'tommo_ac': np.nan,
                                      'tommo_filter': '', 'tommo_rsid': '', 'tommo_match': 'absent'}
                continue
            if len(hits) > 1:
                log(f'    ToMMo: {len(hits)} records for {r.variant_id}, using the first: {hits}')
            rsid, filt, an, ac, af = hits[0]
            direct = (r.OA, r.EA) in raw_alleles(vcf, r.chrom, r.pos)
            rows[r.variant_id] = {'tommo_af': float(af), 'tommo_an': int(an), 'tommo_ac': int(ac),
                                  'tommo_filter': filt, 'tommo_rsid': rsid,
                                  'tommo_match': 'direct' if direct else 'normalised'}
    t = pd.DataFrame.from_dict(rows, orient='index').rename_axis('variant_id').reset_index()
    # cross-check against the project's 60KJPN normalised VCF (INFO/AF)
    af60 = {}
    if crosscheck:
        for _, r in v.iterrows():
            for line in run([TABIX, crosscheck, f'chr{r.chrom}:{r.pos}-{r.pos}']).splitlines():
                f = line.split('\t')
                if f'chr{chrom_plain(f[0])}:{f[1]}:{f[3]}:{f[4]}' == r.variant_id:
                    info = dict(kv.split('=', 1) for kv in f[7].split(';') if '=' in kv)
                    af60[r.variant_id] = float(info.get('AF', 'nan'))
    t['tommo60_af'] = t.variant_id.map(af60)
    both = t.dropna(subset=['tommo_af', 'tommo60_af'])
    r_ = np.corrcoef(both.tommo_af, both.tommo60_af)[0, 1] if len(both) > 2 else np.nan
    log(f'ToMMo: found {t.tommo_af.notna().sum()}/{len(t)} '
        f'(direct {int((t.tommo_match == "direct").sum())}, normalised {int((t.tommo_match == "normalised").sum())}); '
        f'vs 60KJPN on {len(both)}: r = {r_:.4f}, max |delta| = {(both.tommo_af - both.tommo60_af).abs().max():.4f}')
    if len(both) > 2:
        assert r_ > 0.99, f'ToMMo 61KJPN vs 60KJPN r = {r_:.4f}'
    return t


def fetch_kg(v, gt_template, sites_template, eas_bfile, jpt_list, expect_eas, expect_jpt, threads, plink_mb):
    kg = Path('kg')
    kg.mkdir(parents=True, exist_ok=True)
    fam = pd.read_csv(f'{eas_bfile}.fam', sep=r'\s+', header=None, dtype=str)
    eas = sorted(fam[1])
    jpt_all = set(pd.read_csv(jpt_list, sep=r'\s+', header=None, dtype=str)[1])
    jpt = sorted(s for s in eas if s in jpt_all)
    assert len(eas) == expect_eas and len(jpt) == expect_jpt, (len(eas), len(jpt))
    (kg / 'eas.txt').write_text('\n'.join(eas) + '\n')
    (kg / 'jpt.txt').write_text('\n'.join(jpt) + '\n')
    # samples must exist in the NYGC call set
    listed = set(run([BCFTOOLS, 'query', '-l', tpl(gt_template, v.chrom.iloc[0])]).split())
    assert set(eas) <= listed, set(eas) - listed

    rows = {vid: {} for vid in v.variant_id}
    for locus, (chrom, region) in locus_windows(v).items():
        vcf = tpl(gt_template, chrom)
        for tag, sfile in (('eas', kg / 'eas.txt'), ('jpt', kg / 'jpt.txt')):
            recs = bcf_norm_query(vcf, region, '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AN\t%INFO/AC\n',
                                  samples=sfile, keep_info='INFO/AC,INFO/AN')
            by_key = {}
            for c, p, ref, alt, filt, an, ac in recs:
                by_key.setdefault(f'chr{chrom_plain(c)}:{p}:{ref}:{alt}', []).append((filt, int(an), int(ac)))
            for _, r in v[v.locus_id == locus].iterrows():
                hits = by_key.get(r.variant_id, [])
                if not hits:
                    rows[r.variant_id].update({f'kg_{tag}_af': np.nan, f'kg_{tag}_an': np.nan})
                    if tag == 'eas':
                        rows[r.variant_id].update({'kg_filter': '', 'kg_match': 'absent'})
                    continue
                assert len(hits) == 1, (r.variant_id, hits)
                filt, an, ac = hits[0]
                rows[r.variant_id].update({f'kg_{tag}_af': ac / an if an else np.nan, f'kg_{tag}_an': an})
                if tag == 'eas':
                    direct = (r.OA, r.EA) in raw_alleles(vcf, r.chrom, r.pos)
                    rows[r.variant_id].update({'kg_filter': filt, 'kg_match': 'direct' if direct else 'normalised'})
    k = pd.DataFrame.from_dict(rows, orient='index').rename_axis('variant_id').reset_index()

    # cross-check A: plink2 --freq on the project EAS LD panel (same samples)
    ids_file = kg / 'ids.txt'
    ids_file.write_text('\n'.join(sorted(v.variant_id)) + '\n')
    keep = kg / 'jpt.keep'
    keep.write_text(''.join(f'0\t{s}\n' for s in jpt))
    # plink2 sizes its workspace from the node's RAM, not the job's limit, and the
    # 67-M-variant panel needs ~5.7 GB: cap it explicitly (--memory, MB)
    mem = ['--memory', plink_mb]
    run([PLINK2, '--bfile', eas_bfile, '--extract', ids_file, '--freq', '--threads', threads, *mem,
         '--out', kg / 'eas_plink'])
    run([PLINK2, '--bfile', eas_bfile, '--extract', ids_file, '--keep', keep, '--freq', '--threads', threads, *mem,
         '--out', kg / 'jpt_plink'])
    for tag in ('eas', 'jpt'):
        fr = read_tsv(kg / f'{tag}_plink.afreq')
        fr = fr.merge(v[['variant_id', 'EA']], left_on='ID', right_on='variant_id')
        af = np.where(fr.ALT == fr.EA, num(fr.ALT_FREQS), 1 - num(fr.ALT_FREQS))
        k[f'kg_{tag}_af_plink'] = k.variant_id.map(dict(zip(fr.ID, af)))
        both = k.dropna(subset=[f'kg_{tag}_af', f'kg_{tag}_af_plink'])
        # same samples, but the LD panel is another normalisation of the call set
        # and differs by a single allele at a few multi-allelic sites; anything
        # beyond one allele would mean a sample or allele mismatch
        alleles = ((both[f'kg_{tag}_af'] - both[f'kg_{tag}_af_plink']).abs() * both[f'kg_{tag}_an'])
        off = both[alleles > 0.5]
        if not len(both):
            log(f'1KG {tag.upper()}: NYGC route {k[f"kg_{tag}_af"].notna().sum()}/{len(k)}; no variant in both routes')
            continue
        log(f'1KG {tag.upper()}: NYGC route {k[f"kg_{tag}_af"].notna().sum()}/{len(k)}, '
            f'LD-panel plink route {k[f"kg_{tag}_af_plink"].notna().sum()}/{len(k)}, '
            f'{len(off)} of {len(both)} shared variants differ by one allele'
            + (': ' + ', '.join(off.variant_id) if len(off) else ''))
        assert alleles.max() < 1.5, alleles.max()

    # cross-check B: INFO/AF_EAS_unrel of the NYGC sites file at biallelic, directly matching sites
    chk = []
    for _, r in k.merge(v[['variant_id', 'chrom', 'pos', 'OA', 'EA']], on='variant_id').iterrows():
        if r.kg_match != 'direct' or pd.isna(r.kg_eas_af):
            continue
        for line in run([TABIX, tpl(sites_template, r.chrom), f'chr{r.chrom}:{r.pos}-{r.pos}']).splitlines():
            f = line.split('\t')
            alts = f[4].split(',')
            if f[1] != str(r.pos) or f[3] != r.OA or r.EA not in alts or len(alts) != 1:
                continue
            info = dict(kv.split('=', 1) for kv in f[7].split(';') if '=' in kv)
            chk.append(abs(float(info['AF_EAS_unrel']) - r.kg_eas_af))
    if chk:
        log(f'1KG EAS vs INFO/AF_EAS_unrel at {len(chk)} biallelic sites: max |delta| = {max(chk):.2e}')
        assert max(chk) < 1e-3
    else:
        log('1KG EAS vs INFO/AF_EAS_unrel: no biallelic site present (rare variants), check not applicable')
    return k


def main():
    global BCFTOOLS, TABIX, FASTA, PLINK2
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='source', required=True)
    t = sub.add_parser('tommo')
    t.add_argument('--vcf-template', required=True, help='per-chromosome VCF, @@CHROM@@ = 1..22')
    t.add_argument('--crosscheck-vcf', default='', help='60KJPN normalised VCF (optional)')
    k = sub.add_parser('kg')
    k.add_argument('--gt-template', required=True)
    k.add_argument('--sites-template', required=True)
    k.add_argument('--eas-bfile', required=True)
    k.add_argument('--jpt-list', required=True)
    k.add_argument('--expect-eas', type=int, required=True)
    k.add_argument('--expect-jpt', type=int, required=True)
    k.add_argument('--plink2', required=True)
    k.add_argument('--threads', type=int, default=1)
    k.add_argument('--plink-memory', type=int, default=6000, help='plink2 --memory (MB) for the LD-panel cross-check')
    for p in (t, k):
        p.add_argument('--variants', required=True)
        p.add_argument('--out', required=True)
    for p in (t, k):
        p.add_argument('--bcftools', required=True)
        p.add_argument('--tabix', required=True)
        p.add_argument('--fasta', required=True)
    a = ap.parse_args()
    v = unique_variants(a.variants)
    BCFTOOLS, TABIX, FASTA = a.bcftools, a.tabix, a.fasta
    if a.source == 'tommo':
        out = fetch_tommo(v, a.vcf_template, a.crosscheck_vcf)
    else:
        PLINK2 = a.plink2
        out = fetch_kg(v, a.gt_template, a.sites_template, a.eas_bfile, a.jpt_list,
                       a.expect_eas, a.expect_jpt, a.threads, a.plink_memory)
    assert len(out) == len(v) and set(out.variant_id) == set(v.variant_id)
    write_tsv(out, a.out)


if __name__ == '__main__':
    main()
