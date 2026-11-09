#!/usr/bin/env python3
"""Figures of check/platform_depth, six families (one README each, written by catalogue.py):

  01.depth          effective depth by platform and the down-sampling target; local depth at the calling
                    intervals; the read-length audit
  02.downsampling   genotype concordance; the GWAS leads, Tier 1-2 genes and HLA signals across arms
  03.platform       per-platform estimates of the signals with the meta-analyses and leave-one-platform-out;
                    heterogeneity P against the null SNPs
  04.array          the SNP-array evidence per signal, one figure per cohort (11.array)
  05.controls       per-control-source estimates of the signals (12.controls)
  06.evidence       the evidence matrix (14.evidence), z observed against expected for every leave-one-out,
                    and the follow-up counts (13.followup): carriers and allele frequencies by group (one figure
                    per cohort), DQB1:57:V across the replicates

Every number drawn comes from the result tables; the shared project style (analysis/_shared/scripts/
plot_style.py) sets fonts, colours and the cohort ramp.
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import checks as CK  # noqa: E402
import common as K  # noqa: E402
import report_lib as RL  # noqa: E402
sys.path.insert(0, str(K.PROJECT / 'analysis' / '_shared' / 'scripts'))
import plot_style as S  # noqa: E402

PLAT_ORDER = [K.CONTROL_PLATFORM] + K.CASE_PLATFORMS
PLAT_SHORT = RL.PLATFORM_SHORT
COH = dict(zip(K.COHORTS, S.COHORT_RAMP))
COH_SHORT = RL.COHORT_SHORT
STRATA = RL.PLATFORM_ONLY


def save(fig, out, name, meta):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f'{name}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    (out / f'{name}.stats.json').write_text(json.dumps(meta, indent=2, default=float) + '\n')


# ── 01.depth ────────────────────────────────────────────────────────────────
def fig_depth(res, out):
    d = pd.read_csv(res / '01.depth' / 'depth_per_sample.tsv', sep='\t', dtype={'sample_id': str})
    d = d[d.in_full_mainland == 1]
    summ = json.loads((res / '01.depth' / 'depth_summary.json').read_text())
    fa = res / '01.depth' / 'fractions.all.tsv'          # every down-sampled case (T7_FRACTIONS), else the deep cases
    fr = pd.read_csv(fa if fa.exists() else res / '01.depth' / 'fractions.tsv', sep='\t', dtype={'sample_id': str})
    ds_plat = set(fr.platform)
    plat = pd.read_csv(res / '01.depth' / 'platform_depth.tsv', sep='\t')
    plat = plat[plat.cohort == 'full_mainland'].set_index('platform')
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [2.2, 1.3, 1.3]})
    ax = axes[0]
    for i, pf in enumerate(PLAT_ORDER):
        v = d[d.platform == pf].probe_depth.values
        x = i + np.random.default_rng(1).uniform(-0.28, 0.28, len(v))
        ax.scatter(x, v, s=2, color=S.DATA if pf != K.CONTROL_PLATFORM else S.REFERENCE, alpha=0.35, lw=0)
        ax.plot([i - 0.32, i + 0.32], [np.median(v)] * 2, color=S.INK, lw=1.2)
        if pf in ds_plat:
            ax.annotate('', xy=(i + 0.38, summ['target_depth']), xytext=(i + 0.38, np.median(v)),
                        arrowprops=dict(arrowstyle='->', color=S.ACCENT, lw=1.1, shrinkA=0, shrinkB=0))
    ax.axhline(summ['target_depth'], color=S.ACCENT, lw=0.9, ls='--')
    ax.text(-0.45, summ['target_depth'] + 0.6, f"target {summ['target_depth']:.1f}×", color=S.ACCENT,
            va='bottom', ha='left', fontsize=6.5)
    ax.set_xticks(range(len(PLAT_ORDER)), [PLAT_SHORT[p] for p in PLAT_ORDER], fontsize=6.5)
    ax.set_ylabel('Effective depth (×)')
    S.panel_tag(ax, 'a', 'Effective depth per genome')
    ax = axes[1]
    bins = np.linspace(0.3, 1.0, 29)
    groups = [(p, g.fraction.astype(float).values) for p, g in fr.groupby('platform', sort=False)]
    ax.hist([v for _, v in groups], bins=bins, stacked=True, color=[S.SERIES[k] for k in range(len(groups))],
            label=[f'{PLAT_SHORT[p].replace(chr(10), " ")} ({len(v)})' for p, v in groups])
    ax.set_xlabel('Fraction of read pairs kept')
    ax.set_ylabel(f'Down-sampled cases (n = {len(fr)})')
    ax.legend(fontsize=5, frameon=False, loc='upper left')
    S.panel_tag(ax, 'b', 'Fraction kept')
    ax = axes[2]
    xs = np.arange(len(PLAT_ORDER))
    ax.bar(xs - 0.18, [plat.loc[p, 'readlen_median'] for p in PLAT_ORDER], 0.36, color=S.DATA, label='measured')
    ax.bar(xs + 0.18, [plat.loc[p, 'readlen_implied_median'] for p in PLAT_ORDER], 0.36, color=S.NEUTRAL_D,
           label='implied by sheet depth')
    ax.set_xticks(xs, [PLAT_SHORT[p] for p in PLAT_ORDER], fontsize=5.5, rotation=90)
    ax.set_ylabel('Read length (bp)')
    ax.set_ylim(0, 200)
    ax.legend(fontsize=5.5, frameon=False, loc='upper center', ncol=1, bbox_to_anchor=(0.5, 1.0))
    S.panel_tag(ax, 'c', 'Read-length audit')
    fig.tight_layout()
    save(fig, out, 'depth_platforms', {'target_depth': summ['target_depth'], 'downsampled_cases': len(fr),
                                       'by_platform': fr.platform.value_counts().to_dict(),
                                       'platform_median_depth': plat.depth_median.to_dict()})

    loc = pd.read_csv(res / '01.depth' / 'local_depth.tsv', sep='\t')
    fig, ax = plt.subplots(figsize=(7.2, 2.3))
    x = np.arange(len(loc))
    t7f = fr[fr.platform == 'DNBSeq-T7 30x'].fraction.astype(float)
    if len(t7f):                                            # one platform fraction: the local depth scales with it
        loc['median_t7_after'] = loc.median_t7 * float(t7f.iloc[0])
    series = [('median_hiseqx', 'HiSeq X controls', S.REFERENCE), ('median_t7', 'DNBSeq-T7 cases', S.SERIES[2])]
    if len(t7f):
        series.append(('median_t7_after', f'DNBSeq-T7 cases × {float(t7f.iloc[0]):.3f}', S.SERIES[3]))
    series.append(('median_deep_after', 'deep cases after down-sampling', S.ACCENT))
    for col, lab, colr in series:
        ax.plot(x, loc[col], 'o-', ms=2.5, lw=0.8, color=colr, label=lab)
    ax.set_xticks(x, [r.split(':')[0] + ':' + f"{int(r.split(':')[1].split('-')[0]) / 1e6:.2f}M" for r in loc.region],
                  rotation=90, fontsize=6.5)
    ax.set_ylabel('Median effective depth (×)')
    ax.legend(fontsize=6, frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(0.5, 1.08))
    S.panel_tag(ax, 'a', pad=16)              # the legend takes the title slot; the caption names the panel
    fig.tight_layout()
    meta = {'ratio_deep_after_over_controls': loc.deep_after_over_controls.describe().to_dict()}
    if len(t7f):
        meta['ratio_t7_after_over_controls'] = (loc.median_t7_after / loc.median_hiseqx).describe().to_dict()
    save(fig, out, 'depth_local', meta)


# ── 02.downsampling ─────────────────────────────────────────────────────────
def fig_concordance(ds, out):
    c = pd.read_csv(ds / 'concordance' / 'concordance.summary.tsv', sep='\t')
    c = c[c.platform.isna()]
    comps = [('recall_vs_original', 'deep', 'Re-processing\n(baseline vs original)'),
             ('replicate_vs_baseline', 'deep', 'Depth\n(replicates vs baseline)')]
    kinds = [('gwas_signal', 'GWAS signals'), ('gwas_context', 'GWAS context'), ('gene', 'Gene variants')]
    mets = [('concordance', 'Concordance\n(both called)'), ('het_retained', 'Heterozygotes\nretained'),
            ('callrate_test', 'Call rate')]
    fig, axes = plt.subplots(1, len(mets), figsize=(7.2, 2.4), sharey=False)
    for ax, (m, mlab) in zip(axes, mets):
        for j, (comp, role, clab) in enumerate(comps):
            for k, (kind, klab) in enumerate(kinds):
                g = c[(c.comparison == comp) & (c.role == role) & (c.kind == kind) & (c.maf_class != 'all')]
                g = g.groupby('arm')[m].mean() if len(g) else pd.Series(dtype=float)
                if not len(g):
                    continue
                xpos = j * 4 + k
                ax.plot([xpos] * len(g), g.values, 'o', ms=3, color=S.SERIES[k], alpha=0.8,
                        label=klab if (j == 0 and m == mets[0][0]) else None)
        ax.set_xticks([1, 5], [x[2] for x in comps], fontsize=6)
        ax.set_title(mlab, fontsize=7)
    axes[0].legend(fontsize=6, frameon=False, loc='lower left')
    fig.tight_layout()
    save(fig, out, 'concordance', {'rows': len(c)})

    # per platform of the down-sampled cases: the same metrics over every fixed site
    c = pd.read_csv(ds / 'concordance' / 'concordance.summary.tsv', sep='\t')
    c = c[c.platform.notna() & (c.role == 'deep') & (c.kind == 'all') & (c.maf_class == 'all')]
    plats = [p for p in K.CASE_PLATFORMS if p in set(c.platform)]
    fig, axes = plt.subplots(1, len(mets), figsize=(7.2, 2.2), sharey=False)
    meta = {}
    for ax, (m, mlab) in zip(axes, mets):
        for j, (comp, _role, clab) in enumerate(comps):
            for k, pf in enumerate(plats):
                g = c[(c.comparison == comp) & (c.platform == pf)][m]
                if not len(g):
                    continue
                ax.plot([j * (len(plats) + 1) + k] * len(g), g.values, 'o', ms=3, color=S.SERIES[k], alpha=0.8,
                        label=PLAT_SHORT[pf].replace('\n', ' ') if (j == 0 and m == mets[0][0]) else None)
                meta[f'{comp}:{pf}:{m}'] = [float(g.min()), float(g.max())]
        ax.set_xticks([(len(plats) - 1) / 2, len(plats) + 1 + (len(plats) - 1) / 2], [x[2] for x in comps], fontsize=6)
        ax.set_title(mlab, fontsize=7)
    axes[0].legend(fontsize=6, frameon=False, loc='lower left')
    fig.tight_layout()
    save(fig, out, 'concordance_platforms', meta)


def fig_ds_signals(g, out):
    if (g / 'ds.gwas.tsv').exists():
        d = pd.read_csv(g / 'ds.gwas.tsv', sep='\t')
        d = d[(d.imputation == 'best_guess') & (d.role == 'Lead')]
        loci = d.drop_duplicates('unit').sort_values('unit').unit.tolist()
        fig, ax = plt.subplots(figsize=(7.2, 3.0))
        for i, loc in enumerate(loci):
            for j, c in enumerate(K.COHORTS):
                r = d[(d.unit == loc) & (d.cohort == c)]
                if not len(r):
                    continue
                r = r.iloc[0]
                x = i + (j - 1) * 0.25
                ax.plot(x, -np.log10(r.p_original), 'o', ms=3.5, color=COH[c], mfc='white', zorder=3)
                ax.plot(x, -np.log10(r.p_baseline), 's', ms=2.8, color=COH[c], zorder=3)
                ax.plot([x, x], [r.mlog10p_rep_min, r.mlog10p_rep_max], color=COH[c], lw=1.4)
        ax.axhline(5, color=S.REFERENCE, ls=':', lw=0.8)
        ax.axhline(-np.log10(5e-8), color=S.ACCENT, ls='--', lw=0.8)
        ax.set_xticks(range(len(loci)), loci, rotation=60, ha='right', fontsize=6, style='italic')
        ax.set_ylabel(S.NEGLOG10P)
        S.panel_tag(ax, 'a', 'GWAS leads: original (open), baseline (square), replicates (line)')
        fig.tight_layout()
        save(fig, out, 'gwas_leads', {'loci': len(loci)})
    if (g / 'ds.hla_dqb1_57v.tsv').exists():
        h = pd.read_csv(g / 'ds.hla_dqb1_57v.tsv', sep='\t')
        fig, ax = plt.subplots(figsize=(3.4, 2.4))
        for j, c in enumerate(K.COHORTS):
            r = h[h.cohort == c]
            if not len(r):
                continue
            r = r.iloc[0]
            ax.plot(j - 0.15, r.p_over_gec_original, 'o', color=COH[c], mfc='white')
            ax.plot(j, r.p_over_gec_baseline, 's', color=COH[c])
            ax.plot([j + 0.15] * 2, [r.p_over_gec_rep_min, r.p_over_gec_rep_max], color=COH[c], lw=1.5)
        ax.axhline(1, color=S.ACCENT, ls='--', lw=0.8)
        ax.set_yscale('log')
        ax.set_xticks(range(3), [COH_SHORT[c] for c in K.COHORTS])
        ax.set_ylabel('P / GEC threshold')
        S.panel_tag(ax, 'a', 'DQB1:57:V')
        fig.tight_layout()
        save(fig, out, 'hla_dqb1_57v', {})


# ── 03.platform ─────────────────────────────────────────────────────────────
def fig_platform(res, out):
    leads = RL.gwas_leads(RL.read_signals(res))
    rows = [(r.unit, r.variant_id) for r in leads.sort_values('unit').itertuples()] + [(RL.V57, RL.V57)]
    for c in RL.cohorts_with(res, lambda r, c: RL.platform_table(r, c, 'heterogeneity')):
        het = pd.read_csv(RL.platform_table(res, c, 'heterogeneity'), sep='\t')
        st = pd.read_csv(RL.platform_table(res, c, 'strata'), sep='\t', dtype={'variant_id': str})
        fig, ax = plt.subplots(figsize=(7.2, 0.28 * len(rows) + 0.9))
        for i, (lab, vid) in enumerate(rows):
            y = len(rows) - i
            h = het[het.variant_id == vid]
            for k, sname in enumerate(STRATA):
                s_ = st[(st.variant_id == vid) & (st.subset == sname)]
                if len(s_):
                    s_ = s_.iloc[0]
                    ax.plot([s_.or_l95, s_.or_u95], [y + 0.12 * (k - 1.5)] * 2, color=S.SERIES[k], lw=0.8, alpha=0.8)
                    ax.plot(s_['or'], y + 0.12 * (k - 1.5), 'o', ms=2.5, color=S.SERIES[k],
                            label=K.CASE_PLATFORMS[k] if i == 0 else None)
            if len(h):
                h = h.iloc[0]
                ax.plot(np.exp(h.beta_pooled), y + 0.36, 'D', ms=3.5, color=S.INK, label='pooled' if i == 0 else None)
                ax.text(30, y, f"Q P {h.p_Q:.2g}  I² {100 * h.I2:.0f}%", fontsize=5.5, va='center')
        ax.axvline(1, color=S.REFERENCE, lw=0.8)
        ax.set_xscale('log')
        ax.set_xlim(0.05, 120)
        ax.set_yticks([len(rows) - i for i in range(len(rows))], [r[0] for r in rows], fontsize=6, style='italic')
        ax.set_xlabel('Odds ratio (platform cases vs all controls)')
        ax.legend(fontsize=5.5, frameon=False, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.02 + 0.6 / len(rows)))
        S.panel_tag(ax, 'a', COH_SHORT[c] + ' mainland')
        fig.tight_layout()
        save(fig, out, f'forest.{c}', {'markers': len(rows)})


# ── 04.array ────────────────────────────────────────────────────────────────
def fig_array(tables, out):
    f = Path(tables) / 'array.evidence.tsv'
    if not f.exists():
        return
    every = pd.read_csv(f, sep='\t')
    sms = pd.read_csv(Path(tables) / 'array.signals.tsv', sep='\t')
    for c in K.COHORTS:
        fig_array_one(every[every.cohort == c], sms[sms.cohort == c], c, out)


def fig_array_one(ev, sm, cohort, out):
    c = ev[(ev.kind == 'gene') & (ev.evidence != 'none')].dropna(subset=['carriers_wgs']).iloc[::-1]
    genes = len(c) > 0                                   # gene carriers: full mainland only
    fig, axes = plt.subplots(1, 2 if genes else 1, figsize=(7.2 if genes else 3.4, 3.0), squeeze=False,
                             gridspec_kw={'width_ratios': [1.0, 1.25] if genes else [1.0]})
    axes = axes[0]

    ax = axes[0]
    marks = {'typed': 'o', 'imputed': 'D', 'imputed (HIBAG)': '^', 'LD proxy': 's'}
    a = ev[(ev.kind != 'gene') & ev.evidence.isin(list(marks))].dropna(subset=['or_wgs', 'or_array'])
    for e in marks:
        g = a[a.evidence == e]
        if len(g):
            ax.scatter(g.or_wgs, g.or_array, s=11, marker=marks[e], color=S.DATA, lw=0, label=e, zorder=3)
    lim = [0.45, 5.0]
    ax.plot(lim, lim, color=S.REFERENCE, lw=0.8, zorder=1)
    for sig, g in a.groupby('signal'):
        r = g.iloc[0]
        ax.annotate(sig, (r.or_wgs, r.or_array), fontsize=5, xytext=(3, -7), textcoords='offset points',
                    style='normal' if ':' in sig else 'italic')
    ticks = [0.5, 1, 2, 3, 4]
    ax.set_xscale('log')
    ax.set_yscale('log')
    for axis in (ax.xaxis, ax.yaxis):     # after set_*scale, which resets the locators
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        axis.set_major_formatter(matplotlib.ticker.FixedFormatter([f'{t:g}' for t in ticks]))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('OR, WGS genotypes')
    ax.set_ylabel('OR, array data')
    ax.legend(fontsize=5, frameon=False, loc='upper left', handletextpad=0.2)
    S.panel_tag(ax, 'a')                     # a narrow panel: the cohort goes in the figure title (plot_style)
    fig.suptitle(f'{COH_SHORT[cohort]} mainland', fontsize=7, x=0.01, ha='left', fontweight='bold')

    if genes:   # the WGS case carriers of each driver and those the imputed arrays confirm (dosage >= 0.5)
        ax = axes[1]
        y = np.arange(len(c))
        ax.barh(y, c.carriers_wgs, color=S.NEUTRAL, height=0.6, label='WGS case carriers')
        ax.barh(y, c.carriers_found, color=S.DATA, height=0.6, label='confirmed on the imputed arrays')
        ax.set_yticks(y, list(c.marker), fontsize=5)
        ax.grid(axis='y', visible=False)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.set_xlabel('Case carriers')
        ax.legend(fontsize=5, frameon=False, loc='upper center', bbox_to_anchor=(0.35, -0.2), ncol=2)
        S.panel_tag(ax, 'b')
    fig.tight_layout()
    save(fig, out, f'array_evidence.{cohort}', {'association_rows': len(a), 'driver_rows': len(c),
                                      'signals_by_verdict': sm.array_verdict.value_counts().to_dict()})


# ── 05.controls ─────────────────────────────────────────────────────────────
def fig_control_sources(tables, out):
    f = Path(tables) / 'controls.strata.tsv'
    if not f.exists():
        return
    st = pd.read_csv(f, sep='\t', dtype={'variant_id': str})
    het = pd.read_csv(Path(tables) / 'controls.heterogeneity.tsv', sep='\t', dtype={'variant_id': str})
    sm = pd.read_csv(Path(tables) / 'controls.signals.tsv', sep='\t')
    srcs = [s.split('_', 1)[1] for s in sorted(st.subset.unique()) if s.startswith('ctrlonly_')]
    srcs = [s for s in RL.SOURCES if s in srcs]
    label = RL.SOURCE_LABEL
    for c in K.COHORTS:
        rows = sm[sm.cohort == c][['signal', 'marker']].drop_duplicates().values.tolist()
        if not rows:
            continue
        rows = sorted([r for r in rows if r[0] != RL.V57]) + [r for r in rows if r[0] == RL.V57]
        fig, ax = plt.subplots(figsize=(7.2, 0.28 * len(rows) + 0.9))
        for i, (lab, vid) in enumerate(rows):
            y = len(rows) - i
            for k, s in enumerate(srcs):
                x = st[(st.cohort == c) & (st.variant_id == vid) & (st.subset == f'ctrlonly_{s}')]
                if len(x):
                    x = x.iloc[0]
                    yy = y + 0.12 * (k - (len(srcs) - 1) / 2)
                    ax.plot([x.or_l95, x.or_u95], [yy] * 2, color=S.SERIES[k], lw=0.8, alpha=0.8)
                    ax.plot(x['or'], yy, 'o', ms=2.5, color=S.SERIES[k], label=label[s] if i == 0 else None)
            h = het[(het.cohort == c) & (het.variant_id == vid)]
            if len(h):
                h = h.iloc[0]
                ax.plot(np.exp(h.beta_pooled), y + 0.36, 'D', ms=3.5, color=S.INK, label='all controls' if i == 0 else None)
                if np.isfinite(h.get('p_Q', np.nan)):
                    ax.text(30, y, f"Q P {h.p_Q:.2g}  I² {100 * h.I2:.0f}%", fontsize=5.5, va='center')
        ax.axvline(1, color=S.REFERENCE, lw=0.8)
        ax.set_xscale('log')
        ax.set_xlim(0.05, 120)
        ax.set_yticks([len(rows) - i for i in range(len(rows))], [r[0] for r in rows], fontsize=6)
        for t in ax.get_yticklabels():
            t.set_style('normal' if ':' in t.get_text() else 'italic')
        ax.set_xlabel('Odds ratio (all cases vs one control source)')
        ax.legend(fontsize=5.5, frameon=False, ncol=6, loc='upper center', bbox_to_anchor=(0.5, 1.02 + 0.6 / len(rows)))
        S.panel_tag(ax, 'a', COH_SHORT[c] + ' mainland')
        fig.tight_layout()
        save(fig, out, f'sources.{c}', {'markers': len(rows), 'sources': srcs})


# ── 06.evidence ─────────────────────────────────────────────────────────────
EV_LABEL = CK.MATRIX_LABEL                 # the registry's columns, every check but discovery
EV_CHECKS = list(EV_LABEL)
EV_COLOR = {'consistent': S.DATA, 'supported': S.DATA, 'not confirmed': S.NEUTRAL_D,
            'flag': S.ACCENT, 'not testable': S.NEUTRAL, 'not checkable': S.NEUTRAL, 'not applicable': 'white',
            'reported': '#E8EBEE'}
GROUPS_CASE = K.CASE_PLATFORMS
GROUPS_CTRL = RL.SOURCES
SRC_LABEL = RL.SOURCE_LABEL
SET_SHORT = RL.SET_SHORT
# the 13.followup groups drawn, top to bottom: each case platform, all cases, each control source, all controls
FOLLOWUP_ROWS = ([(RL.CASE + p, PLAT_SHORT[p].replace('\n', ' ')) for p in GROUPS_CASE] + [(RL.CASE_ALL, 'all cases')] +
                 [(RL.CONTROL + s, SRC_LABEL[s]) for s in GROUPS_CTRL] + [(RL.CONTROL_ALL, 'all controls')])


def ev_label(signal, kind):
    if kind == 'gene':
        g, st = signal.split(' (')
        return rf'$\it{{{g}}}$ {SET_SHORT[st.rstrip(")")]}'
    return signal


def fig_evidence(tables, out):
    f = Path(tables) / 'evidence.long.tsv'
    if not f.exists():
        return
    ev = pd.read_csv(f, sep='\t')
    rank = {'gwas': 0, 'hla_residue': 1, 'hla_position': 2, 'gene': 3}
    sigs = sorted(ev[['signal', 'kind']].drop_duplicates().values.tolist(), key=lambda s: rank[s[1]])
    n = len(sigs)
    fig, ax = plt.subplots(figsize=(7.2, 0.2 * n + 1.6))
    counts = {}
    for i, (s, k) in enumerate(sigs):
        y = n - 1 - i
        for j, chk in enumerate(EV_CHECKS):
            x = ev[(ev.signal == s) & (ev.check == chk)]
            if not len(x):
                continue
            by = dict(zip(x.cohort, x.reading))
            cells = [(0, 1, x.reading.iloc[0])] if len(x) == 1 else \
                [(m / 3, 1 / 3, by.get(c, 'not applicable')) for m, c in enumerate(K.COHORTS)]
            for x0, w, rd in cells:
                counts[rd] = counts.get(rd, 0) + 1
                ax.add_patch(plt.Rectangle((j + 0.1 + 0.8 * x0, y + 0.14), 0.8 * w - 0.03, 0.72,
                                           facecolor=EV_COLOR.get(rd, S.NEUTRAL), edgecolor=S.NEUTRAL_D, lw=0.3))
    for yb in [sum(1 for s in sigs if rank[s[1]] <= r) for r in (0, 2)]:
        ax.axhline(n - yb, color=S.INK_SOFT, lw=0.6)
    ax.grid(False)
    ax.set_xlim(0, len(EV_CHECKS))
    ax.set_ylim(0, n)
    ax.set_xticks(np.arange(len(EV_CHECKS)) + 0.5, [EV_LABEL[c] for c in EV_CHECKS], fontsize=5.5)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    ax.set_yticks(np.arange(n)[::-1] + 0.5, [ev_label(s, k) for s, k in sigs], fontsize=6)
    for t, (s, k) in zip(ax.get_yticklabels(), sigs):
        t.set_style('italic' if k == 'gwas' else 'normal')
    for sp in ax.spines.values():
        sp.set_visible(False)
    shown = [r for r in EV_COLOR if r in counts]
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=EV_COLOR[r], edgecolor=S.NEUTRAL_D, lw=0.3) for r in shown]
    ax.legend(handles, shown, fontsize=5.5, frameon=False, ncol=len(shown), loc='upper center',
              bbox_to_anchor=(0.5, -0.01), handlelength=1.2, columnspacing=1.0)
    ax.text(0, -0.9 - 0.02 * n, 'Three cells per check: narrow / intermediate / full mainland; the gene sets\' arrays: '
            'one cell (counted in full mainland, which contains the other two).',
            fontsize=5.5, color=S.INK_SOFT, transform=ax.transData)
    fig.tight_layout()
    save(fig, out, 'evidence_matrix', {'signals': n, 'cells': counts})


def fig_z_expected(tables, out):
    f = Path(tables) / 'evidence.loo.tsv'
    if not f.exists():
        return
    lo = pd.read_csv(f, sep='\t').dropna(subset=['z_observed', 'z_expected'])
    lo['x'] = lo.z_expected.abs()
    lo['y'] = lo.z_observed * np.sign(lo.z_expected)
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.9), sharex=True, sharey=True)
    lim = (0, max(7.5, float(lo.x.max()) + 0.3))
    zx, zo = CK.LOO_Z_EXPECTED, CK.LOO_Z_OBSERVED
    flagged = {}
    for r, chk in enumerate(['lopo', 'sources_loo']):
        for col, c in enumerate(K.COHORTS):
            ax = axes[r, col]
            x = lo[(lo.check == chk) & (lo.cohort == c)]
            ax.add_patch(plt.Rectangle((zx, -3), lim[1], 3 + zo, color=S.ACCENT_LT, alpha=0.25, lw=0, zorder=0))
            ax.plot(lim, lim, color=S.REFERENCE, lw=0.7, zorder=1)
            ax.axhline(zo, color=S.REFERENCE, lw=0.5, ls=':')
            ax.axvline(zx, color=S.REFERENCE, lw=0.5, ls=':')
            for kind, mk in [('gwas', 'o'), ('hla_residue', 's'), ('gene', '^')]:
                y = x[x.kind == kind]
                bad = RL.loo_flags(y.z_observed, y.z_expected)      # the evidence-matrix rule
                ax.scatter(y.x[~bad], y.y[~bad], s=9, marker=mk, facecolor='none' if kind == 'gene' else S.DATA,
                           edgecolor=S.DATA, lw=0.6, zorder=3, label={'gwas': 'GWAS lead', 'hla_residue': 'DQB1:57:V',
                                                                      'gene': 'gene (burden)'}[kind])
                ax.scatter(y.x[bad], y.y[bad], s=12, marker=mk, color=S.ACCENT, zorder=4)
                for p in y[bad].itertuples():
                    flagged.setdefault(f'{chk}:{c}', []).append(f'{p.signal} -{p.left_out}')
                    ax.annotate(f'{p.signal}\n−{p.left_out}', (p.x, p.y), fontsize=4.5, xytext=(3, -2),
                                textcoords='offset points', color=S.ACCENT, style='italic' if kind != 'hla_residue' else 'normal')
            ax.set_xlim(lim)
            ax.set_ylim(-3, lim[1])
            if r == 0:
                S.panel_tag(ax, 'abc'[col], COH_SHORT[c])
            else:
                S.panel_tag(ax, 'def'[col], COH_SHORT[c])
                ax.set_xlabel('|z| expected (pooled β / SE without)')
        axes[r, 0].set_ylabel(('Leave one platform out' if chk == 'lopo' else 'Leave one control source out') +
                              '\nz observed (expected direction > 0)')
    axes[0, 2].legend(fontsize=5, frameon=False, loc='upper left')
    fig.tight_layout()
    save(fig, out, 'z_observed_expected', {'points': len(lo), 'flagged': flagged})


def fig_gene_carriers(res, tables, out):
    f = Path(tables) / 'followup.gene_carriers.tsv'
    if not f.exists():
        return
    every = pd.read_csv(f, sep='\t')
    for c in K.COHORTS:
        fig_gene_carriers_one(res, every[every.cohort == c], c, out)


def fig_gene_carriers_one(res, fu, c, out):
    sets = RL.gene_sets(RL.read_signals(res))
    rows = FOLLOWUP_ROWS + [('tommo', 'ToMMo-predicted')]
    ncol = 4
    nrow = int(np.ceil(len(sets) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 3.0 * nrow), sharey=True)
    axes = np.atleast_2d(axes)
    meta = {}
    for ax, (g, st), letter in zip(axes.flat, sets, 'abcdefgh'):
        x = fu[(fu.gene == g) & (fu.stratum == st) & (fu.unit == f'{g} (gene)')].set_index('group')
        for i, (key, _) in enumerate(rows):
            y = len(rows) - 1 - i
            if key == 'tommo':
                v = 100 * x.loc[RL.CASE_ALL, 'tommo_expected_carriers'] / x.loc[RL.CASE_ALL, 'n']
                ax.plot(v, y, 'D', ms=3, color=S.INK)
                meta[f'{g}:{st}:tommo'] = round(v, 3)
                continue
            if key not in x.index:
                continue
            r = x.loc[key]
            col = S.ACCENT if key.startswith(RL.CASE) else S.DATA
            ax.plot(100 * r.carrier_share, y, 'o', ms=3.2 if key.endswith('all') else 2.4, color=col)
            lab = f'{int(r.carriers)}/{int(r.n)}' + (f' (+{int(r.hidden_as_missing)} set missing)'
                                                    if r.hidden_as_missing > 0 else '')
            ax.text(100 * r.carrier_share + 0.25, y, lab, fontsize=4.5, ha='left', va='center', color=S.INK_SOFT)
            meta[f'{g}:{st}:{key}'] = round(100 * r.carrier_share, 3)
        ax.axhline(len(rows) - 5.5, color=S.NEUTRAL, lw=0.6)
        ax.axhline(len(rows) - 11.5, color=S.NEUTRAL, lw=0.6)
        ax.set_xlim(0, 15)
        ax.grid(axis='x', lw=0.4)
        ax.grid(axis='y', visible=False)
        ax.set_xlabel('Carriers (%)')
        S.panel_tag(ax, letter, rf'$\it{{{g}}}$ {SET_SHORT[st]} · {COH_SHORT[c]}', tsize=7)
    for ax in axes.flat[len(sets):]:
        ax.set_visible(False)
    for r in range(nrow):
        axes[r, 0].set_yticks(range(len(rows))[::-1], [lab for _, lab in rows], fontsize=5.5)
    fig.tight_layout()
    save(fig, out, f'gene_carriers.{c}', meta)


def fig_lead_frequencies(tables, out):
    f = Path(tables) / 'followup.gwas_af.tsv'
    if not f.exists():
        return
    every = pd.read_csv(f, sep='\t')
    for c in K.COHORTS:
        fig_lead_frequencies_one(every[every.cohort == c], c, out)


def fig_lead_frequencies_one(fa, c, out):
    leads = list(dict.fromkeys(fa.signal))
    rows = FOLLOWUP_ROWS + [(RL.REF_TOMMO, 'ToMMo 61KJPN'), (RL.REF_JPT, '1KG JPT')]
    ncol = 4
    nrow = int(np.ceil(len(leads) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 1.75 * nrow), sharey=True)
    meta = {}
    for ax, u in zip(axes.flat, leads):
        x = fa[fa.signal == u].set_index('group')
        for i, (key, _) in enumerate(rows):
            if key not in x.index or not np.isfinite(x.loc[key, 'alt_freq']):
                continue
            y = len(rows) - 1 - i
            col = S.ACCENT if key.startswith(RL.CASE) else S.DATA if key.startswith(RL.CONTROL) else S.INK
            ax.plot(x.loc[key, 'alt_freq'], y, 'D' if key.startswith(RL.REFERENCE) else 'o',
                    ms=3 if key.endswith('all') else 2.2, color=col)
        ca, co = x.loc[RL.CASE_ALL, 'alt_freq'], x.loc[RL.CONTROL_ALL, 'alt_freq']
        ax.axvline(co, color=S.DATA, lw=0.5, ls=':')
        ax.axvline(ca, color=S.ACCENT, lw=0.5, ls=':')
        ax.set_title(f'{u} · {COH_SHORT[c]}', fontsize=6.5, style='normal' if ':' in u else 'italic')
        ax.tick_params(labelsize=5)
        meta[u] = {'cases': round(ca, 4), 'controls': round(co, 4)}
    for ax in axes.flat[len(leads):]:
        ax.set_visible(False)
    for r in range(nrow):
        axes[r, 0].set_yticks(range(len(rows))[::-1], [lab for _, lab in rows], fontsize=4.8)
    for ax in axes[-1]:
        ax.set_xlabel('ALT allele frequency', fontsize=5.5)
    fig.tight_layout()
    save(fig, out, f'lead_frequencies.{c}', meta)


def fig_v57_replicates(tables, out):
    f = Path(tables) / 'followup.hla_v57.tsv'
    if not f.exists():
        return
    every = pd.read_csv(f, sep='\t')
    # the deep cases belong to all three cohorts: their typings are the same in each (checked), drawn once
    cols = [k for k in every.columns if k != 'cohort']
    per = [every[every.cohort == c][cols].reset_index(drop=True) for c in K.COHORTS]
    same = all(x.equals(per[-1]) for x in per)
    h = per[-1]
    arms = list(dict.fromkeys(h.arm))
    fig, ax = plt.subplots(figsize=(7.2, 1.9))
    samples = list(dict.fromkeys(h.sample_id))
    for i, s in enumerate(samples):
        x = h[h.sample_id == s].set_index('arm')
        for j, ar in enumerate(arms):
            r = x.loc[ar]
            ch = r.change_vs_baseline
            col = S.ACCENT if ch > 0 else S.DATA if ch < 0 else S.NEUTRAL if np.isfinite(ch) else 'white'
            ax.add_patch(plt.Rectangle((j, len(samples) - 1 - i), 0.94, 0.9, facecolor=col, edgecolor=S.NEUTRAL_D, lw=0.3))
            call = str(r.dqb1_call).replace('DQB1*', '') if isinstance(r.dqb1_call, str) else 'not typed'
            tc = 'white' if ch > 0 else S.INK
            ax.text(j + 0.47, len(samples) - 1 - i + 0.58, call, fontsize=4.2, ha='center', va='center', color=tc)
            v = 'NA' if not np.isfinite(r.v57_dosage) else f'{r.v57_dosage:.0f} V'
            conf = '' if not isinstance(r.confidence, str) else f' · {r.confidence}'
            ax.text(j + 0.47, len(samples) - 1 - i + 0.25, v + conf, fontsize=4.2, ha='center', va='center',
                    color='white' if ch > 0 else S.INK_SOFT)
    ax.grid(False)
    ax.set_xlim(0, len(arms))
    ax.set_ylim(0, len(samples))
    ax.set_xticks(np.arange(len(arms)) + 0.47, arms, fontsize=5.5)
    ax.set_yticks(np.arange(len(samples))[::-1] + 0.45, samples, fontsize=5.5)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.tight_layout()
    save(fig, out, 'v57_replicates', {'deep_cases_changed': len(samples), 'arms': arms,
                                      'identical_in_all_three_cohorts': bool(same)})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results', required=True)
    ap.add_argument('--label', required=True, help='the down-sampling run label')
    ap.add_argument('--out', required=True)
    ap.add_argument('--ds-summary', help='folder of the ds.* tables (default: the published one)')
    ap.add_argument('--array-tables', help='folder of the array.* tables (default: results/11.array)')
    ap.add_argument('--control-tables', help='folder of the controls.* tables (default: results/12.controls)')
    ap.add_argument('--followup-tables', help='folder of the followup.* tables (default: results/13.followup)')
    ap.add_argument('--evidence-tables', help='folder of the evidence.* tables (default: results/14.evidence)')
    a = ap.parse_args()
    S.setup_style('paper')
    res = Path(a.results)
    out = Path(a.out)
    fig_depth(res, out / '01.depth')
    ds = res / '07.downsampling' / a.label
    if (ds / 'concordance').exists():
        fig_concordance(ds, out / '02.downsampling')
    fig_ds_signals(Path(a.ds_summary) if a.ds_summary else ds / 'summary', out / '02.downsampling')
    fig_platform(res, out / '03.platform')
    fig_array(a.array_tables or res / '11.array', out / '04.array')
    fig_control_sources(a.control_tables or res / '12.controls', out / '05.controls')
    ev, fu = a.evidence_tables or res / '14.evidence', a.followup_tables or res / '13.followup'
    fig_evidence(ev, out / '06.evidence')
    fig_z_expected(ev, out / '06.evidence')
    fig_gene_carriers(res, fu, out / '06.evidence')
    fig_lead_frequencies(fu, out / '06.evidence')
    fig_v57_replicates(fu, out / '06.evidence')


if __name__ == '__main__':
    main()
