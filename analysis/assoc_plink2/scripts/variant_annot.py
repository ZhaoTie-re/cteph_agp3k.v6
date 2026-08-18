#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One implementation of every per-variant annotation this component
#           performs, so the lead table, the cross-cohort table and the
#           credible-set table cannot drift apart:
#             rsID     a tabix-indexed VCF carrying rsIDs, keyed chr:pos:REF:ALT
#             gene     EnsDb.Hsapiens.v86 via gene_annotate.R (nearest/overlapping)
#             snpEff   the platform's PRE-COMPUTED snpEff index — the same
#                      resource analysis/assoc_rvtest annotates from
#             genotype per-group counts, EAF, missing rate and HWE from two
#                      plink2 calls (--keep cases / --keep controls)
#           `annotate_variants()` assembles all four into the one row shape that
#           every variant-level table in this component uses. Everything is a
#           tabix / plink2 / R call against resources already in the project;
#           nothing here runs snpEff or queries a network.
# Component: assoc_plink2
# ---------------------------------------------------------------------------
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# snpEff index layout: one tabix-indexed TSV per chromosome, 35 columns.
SNPEFF_TEMPLATE = 'all.VQSR3.chr{chrom}.vcf_out.tsv.2.gz'
SNPEFF_COLS = {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'effect': 8, 'impact': 9,
               'gene': 10, 'geneid': 11, 'feature': 12, 'featureid': 13,
               'biotype': 14, 'hgvs_c': 16, 'hgvs_p': 17}
# snpEff's own ordering; anything unrecognised sorts last.
IMPACT_RANK = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}

SNPEFF_FIELDS = ['snpEff_Effect', 'snpEff_Impact', 'snpEff_Gene', 'snpEff_GeneID',
                 'snpEff_Feature', 'snpEff_FeatureID', 'snpEff_Biotype',
                 'snpEff_HGVS_c', 'snpEff_HGVS_p', 'snpEff_n_annotations',
                 'snpEff_all_effects']
SNPEFF_MISSING = {k: ('.' if k != 'snpEff_n_annotations' else 0) for k in SNPEFF_FIELDS}

# The one variant-level field set. Every table that describes a variant in a
# cohort carries exactly these, in this order, after its own identifying prefix
# (cohort / peak_id / tier, or cohort / model / hit_id). Defined once so the lead
# table, the per-model hit table and the cross-cohort table are directly
# stackable and a reader learns the columns only once.
# OR, L95 and U95 are THREE NUMERIC COLUMNS, never one formatted string. A packed
# '1.591 (1.350–1.874)' cannot be used without parsing it back out, which the
# comparison figure really did (split('(')[1].rstrip(')').split('–')), and it made
# an en-dash a load-bearing character. Format for display at the point of display.
STAT_FIELDS = [
    'rsID', 'Gene', 'Gene_Biotype', 'Gene_Distance_bp',
    'EA', 'OA', 'Beta', 'SE', 'OR', 'L95', 'U95', 'P',
    'Case_Genotype_Distribution', 'Case_EAF', 'Case_Missing_Rate', 'Case_HWE_P',
    'Control_Genotype_Distribution', 'Control_EAF', 'Control_Missing_Rate', 'Control_HWE_P',
    'A1_FREQ', 'MAF', 'OBS_CT', 'N_case', 'N_ctrl',
] + SNPEFF_FIELDS


def or_ci_text(orv, lo, hi):
    """Human-readable 'OR (L95-U95)' for a caption or a printed table.

    The ONLY place the packed form is allowed to exist. Stored tables keep the
    three numeric columns.
    """
    try:
        orv, lo, hi = float(orv), float(lo), float(hi)
    except (TypeError, ValueError):
        return '.'
    if not (np.isfinite(orv) and np.isfinite(lo) and np.isfinite(hi)):
        return '.'
    return f'{orv:.3f} ({lo:.3f}\u2013{hi:.3f})'


def _plain(chrom):
    return str(chrom).replace('chr', '')


def lookup_rsids(tabix, vcf, variants, timeout=120):
    """rsID per variant from the configured rsID VCF.

    `variants` is an iterable of (variant_id, chrom, pos). The VCF is normalised
    and its ID column is keyed exactly like our chr:pos:REF:ALT variant IDs, so
    the match is on the reconstructed key rather than on position alone —
    position alone would assign a multi-allelic site's rsID to the wrong allele.
    """
    out = {}
    if not vcf or not Path(vcf).exists():
        return out
    for vid, chrom, pos in variants:
        region = f'chr{_plain(chrom)}:{int(pos)}-{int(pos)}'
        try:
            res = subprocess.run([tabix, vcf, region], capture_output=True,
                                 text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            continue
        for line in res.stdout.splitlines():
            f = line.split('\t')
            if len(f) < 5:
                continue
            if f'chr{_plain(f[0])}:{f[1]}:{f[3]}:{f[4]}' == vid and f[2] not in ('.', ''):
                out[vid] = f[2]
                break
    return out


def lookup_genes(rscript, gene_script, positions, work_dir):
    """Nearest/overlapping gene per position, via gene_annotate.R (EnsDb v86).

    `positions` is an iterable of (chrom, pos). Returns {(chrom, pos): Series}.
    """
    out = {}
    if not (rscript and gene_script and Path(gene_script).exists()):
        return out
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    pos_file, out_file = work / 'annot_positions.tsv', work / 'annot_genes.tsv'
    pd.DataFrame(sorted({(str(c), int(p)) for c, p in positions}),
                 columns=['chrom', 'pos']).to_csv(pos_file, sep='\t', index=False)
    res = subprocess.run([rscript, gene_script, '--positions', str(pos_file),
                          '--out', str(out_file)], capture_output=True, text=True)
    if res.returncode != 0 or not out_file.exists():
        print((res.stdout + res.stderr)[-2000:], file=sys.stderr)
        return out
    g = pd.read_csv(out_file, sep='\t', dtype={'chrom': str})
    return {(str(r['chrom']), int(r['pos'])): r for _, r in g.iterrows()}


def lookup_snpeff(tabix, index_dir, variants, timeout=120):
    """Functional consequence per variant from the pre-computed snpEff index.

    A variant carries ONE ROW PER TRANSCRIPT, so the rows are collapsed to the
    most severe by IMPACT, preferring a protein_coding feature on ties. The
    discarded rows are not thrown away silently: `snpEff_n_annotations` records
    how many there were and `snpEff_all_effects` lists the distinct effects, so a
    reader can see that a MODIFIER call had, say, a MODERATE sibling transcript.
    """
    out = {}
    if not index_dir or not Path(index_dir).exists():
        return out
    for vid, chrom, pos in variants:
        c = _plain(chrom)
        path = Path(index_dir) / SNPEFF_TEMPLATE.format(chrom=c)
        if not path.exists():
            continue
        try:
            res = subprocess.run([tabix, str(path), f'chr{c}:{int(pos)}-{int(pos)}'],
                                 capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            continue
        rows = []
        for line in res.stdout.splitlines():
            f = line.split('\t')
            if len(f) <= max(SNPEFF_COLS.values()):
                continue
            g = {k: f[i] for k, i in SNPEFF_COLS.items()}
            if f'chr{_plain(g["chrom"])}:{g["pos"]}:{g["ref"]}:{g["alt"]}' == vid:
                rows.append(g)
        if not rows:
            continue
        rows.sort(key=lambda g: (IMPACT_RANK.get(g['impact'], 9),
                                 0 if g['biotype'] == 'protein_coding' else 1))
        best = rows[0]
        effects = sorted({g['effect'] for g in rows if g['effect'] not in ('', '.')})
        out[vid] = {
            'snpEff_Effect': best['effect'] or '.', 'snpEff_Impact': best['impact'] or '.',
            'snpEff_Gene': best['gene'] or '.', 'snpEff_GeneID': best['geneid'] or '.',
            'snpEff_Feature': best['feature'] or '.', 'snpEff_FeatureID': best['featureid'] or '.',
            'snpEff_Biotype': best['biotype'] or '.', 'snpEff_HGVS_c': best['hgvs_c'] or '.',
            'snpEff_HGVS_p': best['hgvs_p'] or '.', 'snpEff_n_annotations': len(rows),
            'snpEff_all_effects': ','.join(effects) or '.',
        }
    return out


# ── genotype-level statistics, one group at a time ───────────────────────────
def _sh(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        print(' '.join(str(c) for c in cmd), file=sys.stderr)
        print((r.stdout + r.stderr)[-3000:], file=sys.stderr)
    return r.returncode == 0


def group_stats(plink2, bfile, fam, work, group, pheno_code, ids, threads=1):
    """Genotype counts, EAF, missing rate and HWE for one phenotype group.

    plink2 computes both reports over whatever `--keep` leaves in the sample set,
    so restricting to cases gives case-only statistics; no phenotype needs to be
    loaded for either report.

    Returns ({variant_id: stats}, n_in_group).
    """
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    keep = work / f'{group}.keep'
    sel = fam[fam['PHENO'] == pheno_code][['FID', 'IID']]
    sel.to_csv(keep, sep='\t', header=False, index=False)
    ex = work / f'{group}.extract'
    ex.write_text('\n'.join(map(str, ids)) + '\n')

    pfx = work / group
    ok = _sh([plink2, '--bfile', bfile, '--keep', keep, '--extract', ex,
              '--geno-counts', '--hardy', '--threads', threads, '--out', pfx])
    n_group = len(sel)
    if not ok:
        return {}, n_group

    gc = pd.read_csv(f'{pfx}.gcount', sep='\t', dtype={'ID': str}) \
        if Path(f'{pfx}.gcount').exists() else pd.DataFrame()
    hw = pd.read_csv(f'{pfx}.hardy', sep='\t', dtype={'ID': str}) \
        if Path(f'{pfx}.hardy').exists() else pd.DataFrame()
    hwe = dict(zip(hw['ID'], hw['P'])) if len(hw) else {}

    out = {}
    for _, r in gc.iterrows():
        hom_ref = float(r.get('HOM_REF_CT', np.nan))
        het = float(r.get('HET_REF_ALT_CTS', np.nan))
        hom_alt = float(r.get('TWO_ALT_GENO_CTS', np.nan))
        miss = float(r.get('MISSING_CT', 0) or 0)
        called = hom_ref + het + hom_alt
        out[r['ID']] = {
            'ref': r.get('REF', ''), 'alt': r.get('ALT', ''),
            'geno': f'{int(hom_ref)}/{int(het)}/{int(hom_alt)}',
            'hom_ref': hom_ref, 'het': het, 'hom_alt': hom_alt, 'called': called,
            'missing_rate': (miss / (called + miss)) if (called + miss) > 0 else np.nan,
            'hwe_p': hwe.get(r['ID'], np.nan),
        }
    return out, n_group


def annotate_variants(records, *, bfile, plink2, tabix, work, threads=1,
                      rsid_vcf=None, rscript=None, gene_script=None, snpeff_index=None):
    """The complete variant-level row for every record, in `STAT_FIELDS` order.

    `records` is an iterable of mappings with the association statistics already
    attached: variant_id, chrom, pos, a1, ref, alt, or, se, l95, u95, p, a1_freq,
    obs_ct. Everything else — rsID, gene, snpEff, and the two per-group genotype
    reports — is looked up here, once, for the whole set.

    Returns a DataFrame indexed positionally over `records`, with columns
    ['chrom', 'pos', 'variant_id'] + STAT_FIELDS. Callers prepend their own
    identifying columns (cohort/peak_id/tier, or cohort/model/hit_id).
    """
    recs = [dict(r) for r in records]
    if not recs:
        return pd.DataFrame(columns=['chrom', 'pos', 'variant_id'] + STAT_FIELDS)

    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    fam = pd.read_csv(f'{bfile}.fam', sep=r'\s+', header=None,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'],
                      dtype={'FID': str, 'IID': str})
    ids = sorted({str(r['variant_id']) for r in recs})
    case, n_case = group_stats(plink2, bfile, fam, work, 'case', 2, ids, threads)
    ctrl, n_ctrl = group_stats(plink2, bfile, fam, work, 'control', 1, ids, threads)

    keys = sorted({(str(r['variant_id']), str(r['chrom']), int(r['pos'])) for r in recs})
    rsids = lookup_rsids(tabix, rsid_vcf, keys)
    snpeff = lookup_snpeff(tabix, snpeff_index, keys)
    genes = lookup_genes(rscript, gene_script, [(c, p) for _v, c, p in keys], work)

    def grp(d, vid, ea):
        """EAF is reported for the EFFECT allele, so the ALT-based count must be
        flipped when plink2 tested REF."""
        g = d.get(vid)
        if not g or not g['called']:
            return '.', np.nan, np.nan, np.nan
        alt_ac = 2 * g['hom_alt'] + g['het']
        ref_ac = 2 * g['hom_ref'] + g['het']
        eaf = (alt_ac if ea == g['alt'] else ref_ac) / (2 * g['called'])
        return g['geno'], round(eaf, 5), round(g['missing_rate'], 5), g['hwe_p']

    rows = []
    for r in recs:
        vid = str(r['variant_id'])
        ea = str(r.get('a1', '') or '')
        ref, alt = str(r.get('ref', '') or ''), str(r.get('alt', '') or '')
        oa = ref if ea == alt else alt
        orv = float(r['or']) if pd.notna(r.get('or')) else np.nan
        se = r.get('se', np.nan)
        lo, hi = r.get('l95', np.nan), r.get('u95', np.nan)
        af = float(r['a1_freq']) if pd.notna(r.get('a1_freq')) else np.nan
        c_geno, c_eaf, c_miss, c_hwe = grp(case, vid, ea)
        k_geno, k_eaf, k_miss, k_hwe = grp(ctrl, vid, ea)
        gi = genes.get((str(r['chrom']), int(r['pos'])))
        row = {
            'chrom': r['chrom'], 'pos': int(r['pos']), 'variant_id': vid,
            'rsID': rsids.get(vid, '.'),
            'Gene': (gi['Gene'] if gi is not None else '.'),
            'Gene_Biotype': (gi['Gene_Biotype'] if gi is not None else '.'),
            'Gene_Distance_bp': (gi['Gene_Distance_bp'] if gi is not None else np.nan),
            'EA': ea, 'OA': oa,
            'Beta': round(float(np.log(orv)), 5) if np.isfinite(orv) and orv > 0 else np.nan,
            'SE': round(float(se), 5) if pd.notna(se) else np.nan,
            'OR': round(float(orv), 5) if np.isfinite(orv) else np.nan,
            'L95': round(float(lo), 5) if pd.notna(lo) else np.nan,
            'U95': round(float(hi), 5) if pd.notna(hi) else np.nan,
            'P': r.get('p', np.nan),
            'Case_Genotype_Distribution': c_geno, 'Case_EAF': c_eaf,
            'Case_Missing_Rate': c_miss, 'Case_HWE_P': c_hwe,
            'Control_Genotype_Distribution': k_geno, 'Control_EAF': k_eaf,
            'Control_Missing_Rate': k_miss, 'Control_HWE_P': k_hwe,
            'A1_FREQ': af, 'MAF': round(min(af, 1 - af), 5) if np.isfinite(af) else np.nan,
            'OBS_CT': r.get('obs_ct', np.nan),
            'N_case': n_case, 'N_ctrl': n_ctrl,
        }
        # snpEff names the transcript the variant actually falls in; `Gene` above
        # is the EnsDb nearest/overlapping gene. Different questions, both kept.
        row.update(snpeff.get(vid, SNPEFF_MISSING))
        rows.append(row)

    for f in list(work.glob('*.gcount')) + list(work.glob('*.hardy')) + list(work.glob('*.log')):
        f.unlink()
    return pd.DataFrame(rows, columns=['chrom', 'pos', 'variant_id'] + STAT_FIELDS)
