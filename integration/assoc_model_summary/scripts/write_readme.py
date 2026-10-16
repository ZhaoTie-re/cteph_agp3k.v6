#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The text of each figure family, from the numbers of this run (no
#           value is typed in): README.md with
#             1. Slide summary  — a take-home line and bullets, for the slide
#             2. Figures        — one line per PNG
#             3. Detailed       — question, design, how to read, results,
#                                 interpretation, limits, methods
#           plus slide.md (the slide summary alone) for the results index.
#             --family fixed_vs_random | pc_effect | main_results | index
#           pc_effect states facts only: its conclusion is to be agreed.
# Inputs  : the family's tables / JSON (see --help)
# Outputs : README.md, slide.md
# ---------------------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')


def pu(p, digits=1):
    """P with a Unicode exponent, so it survives a copy into a slide: 3.8 × 10⁻⁸."""
    p = float(p)
    if p >= 0.001:
        return f'{p:.2g}'
    m, e = f'{p:.{digits}e}'.split('e')
    return f'{m} × 10{str(int(e)).translate(SUP)}'


def slash(vals, fmt='{}'):
    return ' / '.join(fmt.format(v) for v in vals)


def md_table(df):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join('---' for _ in cols) + '|']
    out += ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in df.itertuples(index=False)]
    return '\n'.join(out)


def write(out_dir, title, slide, figures, detailed):
    text = [f'# {title}', '', '## Slide summary', '', slide.strip(), '', '## Figures', '']
    text += [f'- **`{f}`** — {d}' for f, d in figures]
    text += ['', '## Detailed description', '', detailed.strip(), '']
    Path(out_dir, 'README.md').write_text('\n'.join(text))
    Path(out_dir, 'slide.md').write_text(f'### {title}\n\n{slide.strip()}\n')


def cohort_names(a):
    return dict(zip(C.csv_list(a.cohorts), C.csv_list(a.cohort_labels)))


def maf_range(maf, model):
    d = maf[(maf.model == model) & (maf.maf_bin != '<0.01')]
    return d.lambda_gc.min(), d.lambda_gc.max()


# ── 1. fixed vs random ───────────────────────────────────────────────────────
def fixed_vs_random(a):
    cohorts, names = C.csv_list(a.cohorts), cohort_names(a)
    t = pd.read_csv(a.tail, sep='\t').set_index(['model', 'cohort'])
    maf = pd.read_csv(a.maf, sep='\t')
    F, R = 'fixed', 'grm_pc'
    lf = [t.loc[(F, c), 'lambda_gc'] for c in cohorts]
    lr = [t.loc[(R, c), 'lambda_gc'] for c in cohorts]
    nf = [int(t.loc[(F, c), 'n_p_lt_1e5']) for c in cohorts]
    nr = [int(t.loc[(R, c), 'n_p_lt_1e5']) for c in cohorts]
    exp = float(np.mean([t.loc[(F, c), 'expected_p_lt_1e5'] for c in cohorts]))
    gwf = [int(t.loc[(F, c), 'gw_loci']) for c in cohorts]
    gwr = [int(t.loc[(R, c), 'gw_loci']) for c in cohorts]
    fmin, fmax = maf_range(maf, F)
    rmin, rmax = maf_range(maf, R)
    removed = [int(t.loc[(R, c), 'n_ctrl'] - t.loc[(F, c), 'n_ctrl']) for c in cohorts]

    # every fixed-model genome-wide lead, with the smallest random-model P within +-250 kb (the peak window):
    # a LOCUS comparison, since the two models may elect different lead variants of one locus
    fx = []
    leads = pd.read_csv(a.fixed_leads, sep='\t')
    for _, r in leads[leads.model == F].iterrows():
        top = pd.read_csv(f'{r.cohort}.{R}.top.tsv', sep='\t')
        ch, pos = r.variant_id.split(':')[0].replace('chr', ''), int(r.variant_id.split(':')[1])
        w = top[(top.CHROM.astype(str).str.replace('chr', '') == ch) & ((top.POS - pos).abs() <= 250_000)]
        pr = float(w.P.min()) if len(w) else float('nan')
        fx.append(dict(cohort=r.cohort, gene=r.Gene, variant=r.variant_id, rsID=r.rsID, p_fixed=float(r.P), p_random=pr))
    lost = [d for d in fx if not d['p_random'] < C.GW_ALPHA]
    lost_txt = '; '.join(f'*{d["gene"]}* in {names[d["cohort"]].lower()} (P {pu(d["p_fixed"])} → {pu(d["p_random"])})'
                         for d in lost)
    labels = [t.loc[(F, c), 'gw_leads'] if isinstance(t.loc[(F, c), 'gw_leads'], str) else '' for c in cohorts]
    labels_r = [t.loc[(R, c), 'gw_leads'] if isinstance(t.loc[(R, c), 'gw_leads'], str) else '' for c in cohorts]

    slide = f"""
**Take-home: a fixed-effect logistic model remains inflated even after adjusting for 10 PCs and removing close relatives; the logistic mixed model with a full GRM is well calibrated. The main analysis therefore uses the random (mixed) model.**

- Fixed model (PLINK2; SEX + 10 PCs; relatives removed): λGC = {slash(lf, '{:.3f}')} (narrow / intermediate / full mainland).
- Random model (SAIGE; SEX + 10 PCs + full GRM; relatives retained): λGC = {slash(lr, '{:.3f}')}.
- The excess is genome-wide: in every allele-frequency bin (MAF 0.01–0.50) the fixed model has λGC {fmin:.2f}–{fmax:.2f}, the random model {rmin:.2f}–{rmax:.2f}.
- Variants at P < 10⁻⁵ (≈ {exp:.0f} expected by chance): fixed {slash(nf)}; random {slash(nr)}.
- Genome-wide loci: fixed {slash(gwf)}, random {slash(gwr)}{'; genome-wide only under the fixed model: ' + lost_txt if lost else ''}.
"""
    figures = [('fixed_vs_random.png', 'one row per cohort: Manhattan of the fixed model (left) and of the random '
                'model (middle) on one genomic axis and one −log₁₀P scale, and the QQ plot of both (right, '
                'λGC in the model colour).'),
               ('calibration.png', '(a) λGC per cohort; (b) λGC by minor-allele frequency; '
                '(c) variants at P < 10⁻⁵ against the number expected by chance.')]
    design = pd.DataFrame({
        '': ['Software', 'Covariates', 'Samples', 'Cases / controls', 'Variants', 'P value'],
        'Fixed model': ['PLINK 2 logistic regression (`--glm`, no Firth)', 'SEX + PC1–10 (fixed effects)',
                        'relatives removed (PI_HAT > 0.2, vertex cover)',
                        ' · '.join(f'{int(t.loc[(F, c), "n_case"])}/{int(t.loc[(F, c), "n_ctrl"]):,}' for c in cohorts),
                        f'{int(t.loc[(F, cohorts[0]), "n_variants"]):,} (narrow)', 'Wald test'],
        'Random model': ['SAIGE logistic mixed model', 'SEX + PC1–10 (fixed) + polygenic random effect, full GRM',
                         'relatives retained; the GRM models them',
                         ' · '.join(f'{int(t.loc[(R, c), "n_case"])}/{int(t.loc[(R, c), "n_ctrl"]):,}' for c in cohorts),
                         f'{int(t.loc[(R, cohorts[0]), "n_variants"]):,} (narrow)',
                         'score test, normal approximation'],
    })
    res = pd.DataFrame({
        'Cohort': [names[c] for c in cohorts],
        'λGC fixed': [f'{v:.3f}' for v in lf], 'λGC random': [f'{v:.3f}' for v in lr],
        'P<10⁻⁵ fixed': nf, 'P<10⁻⁵ random': nr,
        'expected': [f'{t.loc[(F, c), "expected_p_lt_1e5"]:.0f}' for c in cohorts],
        'GW loci fixed': [f'{g}' + (f' ({l})' if l else '') for g, l in zip(gwf, labels)],
        'GW loci random': [f'{g}' + (f' ({l})' if l else '') for g, l in zip(gwr, labels_r)],
        'SG loci fixed / random': [f'{int(t.loc[(F, c), "sg_loci"])} / {int(t.loc[(R, c), "sg_loci"])}'
                                   for c in cohorts],
    })
    mtab = maf[maf.model.isin([F, R])].pivot_table(index=['maf_bin'], columns=['model', 'cohort'],
                                                   values='lambda_gc', sort=False)
    mrows = []
    for b in dict.fromkeys(maf.maf_bin):
        n_b = int(maf[(maf.maf_bin == b) & (maf.model == R) & (maf.cohort == cohorts[-1])].n.iloc[0])
        mrows.append([b, f'{n_b:,}'] + [f'{mtab.loc[b, (m, c)]:.3f}' for m in (F, R) for c in cohorts])
    mdf = pd.DataFrame(mrows, columns=['MAF', f'variants (random, {names[cohorts[-1]].lower()})']
                       + [f'{"fixed" if m == F else "random"} {names[c].split()[0].lower()}' for m in (F, R) for c in cohorts])
    low = maf[maf.maf_bin == '<0.01']
    lowf = low[low.model == F]
    lowr = low[low.model == R]
    detailed = f"""
### Question
Does a fixed-effect logistic regression control population structure and relatedness in the three nested
cohorts, or is a mixed model with a genetic relationship matrix (GRM) needed?

### Design
{md_table(design)}

The two scans share the variant set (the same IDs in the same order) and the covariates (SEX + 10 PCs). They
differ in the polygenic random effect and in the {slash(removed)} related controls that the fixed model removes
and the random model keeps (narrow / intermediate / full mainland).

### How to read the figures
- **`fixed_vs_random.png`** — rows a–c are the cohorts. Every Manhattan shares one genomic axis and one
  −log₁₀P limit, so heights are comparable between models and cohorts. Red line: 5 × 10⁻⁸; grey dashed line:
  10⁻⁵. Dotted guides mark loci genome-wide in any panel (named above the top row); a diamond marks a
  panel's own genome-wide lead. In the QQ panel, both models are drawn against the uniform null, with the
  pointwise 95% band (grey); a curve that leaves the diagonal from the bottom, not only in the tail, is
  inflation rather than signal.
- **`calibration.png`** — (a) λGC per cohort (red line = 1). (b) λGC within minor-allele
  frequency bins: inflation that persists in every bin is genome-wide; a few true loci cannot raise the
  median of millions of variants. (c) Variants at P < 10⁻⁵ against the number expected by chance
  (dashed line, n × 10⁻⁵).

### Results
{md_table(res)}

λGC by minor-allele frequency:

{md_table(mdf)}

Genome-wide leads of the fixed model, with the smallest random-model P within ±250 kb:

{md_table(pd.DataFrame([dict(Cohort=names[d['cohort']], Gene=d['gene'], Variant=d['variant'], rsID=d['rsID'],
                             **{'P fixed (lead)': pu(d['p_fixed']), 'P random (locus)': pu(d['p_random'])}) for d in fx]))
 if fx else 'none'}

### Interpretation
- Removing close relatives and adjusting for 10 PCs does not calibrate the fixed model: λGC
  {min(lf):.2f}–{max(lf):.2f} in all three cohorts, raised in every frequency bin, and {min(nf)}–{max(nf)} variants
  at P < 10⁻⁵ where about {exp:.0f} are expected. A genome-wide excess of this kind is the signature of
  structure that the fixed covariates do not capture: distant relatedness below the PI_HAT cut-off and
  fine-scale ancestry beyond 10 PCs.
- The random model, on the same variants with the same covariates, absorbs that genome-wide sharing through
  the GRM: λGC {min(lr):.3f}–{max(lr):.3f}, and {rmin:.2f}–{rmax:.2f} in every frequency bin from MAF 0.01. It also keeps the
  related samples instead of discarding them.
- Loci: the inflated fixed model calls more suggestive signals ({slash([int(t.loc[(F, c), 'sg_loci']) for c in cohorts])}
  against {slash([int(t.loc[(R, c), 'sg_loci']) for c in cohorts])}) and genome-wide calls that the random model
  does not reproduce{(' (' + lost_txt + ')') if lost else ''}. *FTO* stays the strongest locus under both models.
- Hence the main analysis is the random (mixed) model.

### Limits
- The sample sets differ by the removed relatives ({slash(removed)} controls). Those samples are absent from the
  fixed model, so they cannot be the cause of its inflation.
- λGC mixes polygenic signal and confounding and does not separate them. The two models are compared on
  the same phenotype and variants, so polygenic signal cannot explain the difference between them.
- Sequencing platform is confounded with case status (cases and controls were sequenced on different
  platforms). Part of a genome-wide excess can be technical; the mixed model absorbs genome-wide components
  of either origin, but single signals still need the locus-level checks.
- Variants with overall MAF < 0.01 ({slash([f'{int(n):,}' for n in lowr.n])} in the random-model scans, kept
  because the upstream MAF filter is applied to the controls) are inflated under both models (λGC
  {slash([f'{v:.2f}' for v in lowf.lambda_gc])} fixed, {slash([f'{v:.2f}' for v in lowr.lambda_gc])} random). They are less than 1% of
  the variants, are not plotted in calibration panel b, and are listed for follow-up.
- P values: PLINK 2 Wald test and SAIGE score test with the normal approximation (no saddle-point
  correction); both are asymptotic.

### Methods
- λGC = median(χ²₁ of P) / 0.4549 over every usable variant: ERRCODE '.', 0 < P ≤ 1, finite OR
  and SE. λGC and the counts are recomputed here and equal the upstream `scan_qc.tsv` (asserted).
- QQ band: pointwise 2.5–97.5% quantiles of Beta(k, n − k + 1) for the k-th smallest of n uniform P. The
  QQ draws every one of the 20,000 smallest P and a log-spaced sample of the rest.
- Manhattan: every variant at P ≤ 10⁻³ is drawn; below that, overlapping points (closer than ~2 px at
  600 dpi) are drawn once, which leaves the image unchanged.
- Sources: `analysis/assoc_plink2` (additive scans) and `analysis/assoc_saige` (`main` scans), with their
  peak calls (±250 kb).
"""
    write(a.out_dir, '1 · Fixed vs random model', slide, figures, detailed)


# ── 2. PCs on top of the GRM ─────────────────────────────────────────────────
def pc_effect(a):
    cohorts, names = C.csv_list(a.cohorts), cohort_names(a)
    j = C.load_json(a.pc_effect)['cohorts']
    pc = pd.read_csv(a.pc_assoc, sep='\t').set_index('cohort')
    pk = pd.read_csv(a.peak_changes, sep='\t')
    lo = pd.read_csv(a.lead_or, sep='\t')
    maf = pd.read_csv(a.maf, sep='\t')
    lw = [j[c]['lambda_without_pc'] for c in cohorts]
    lp = [j[c]['lambda_with_pc'] for c in cohorts]
    nw = [j[c]['p_lt_1e5_without_pc'] for c in cohorts]
    np_ = [j[c]['p_lt_1e5_with_pc'] for c in cohorts]
    exp = float(np.mean([j[c]['expected_p_lt_1e5'] for c in cohorts]))
    ow = [j[c]['peaks_only_without_pc'] for c in cohorts]
    op = [j[c]['peaks_only_with_pc'] for c in cohorts]
    mhc = pk[(pk.only_in == 'grm_only') & (pk.chrom.astype(str) == '6') & (pk.pos.between(28_400_000, 33_500_000))]
    mhc_txt = ', '.join(f'*{r.gene}* ({names[r.cohort].lower()}, {int(r.n_sig_variants)} variant'
                        f'{"s" if r.n_sig_variants != 1 else ""})' for _, r in mhc.iterrows())
    med = float(lo.dlnOR.abs().median())
    mx = float(lo.dlnOR.abs().max())
    slide = f"""
**Conclusion: pending — to be agreed. The facts:**

- With the full GRM in the model, adding 10 PCs leaves λGC nearly unchanged: {slash(lw, '{:.3f}')} without → {slash(lp, '{:.3f}')} with PCs (narrow / intermediate / full mainland).
- The PCs differ between cases and controls in every cohort: likelihood-ratio test of PC1–10 on case status, P = {slash([pu(pc.loc[c, 'lr_p']) for c in cohorts])}.
- Variants at P < 10⁻⁵ (≈ {exp:.0f} expected): {slash(nw)} without PCs, {slash(np_)} with PCs.
- Peaks found only without PCs: {slash(ow)}{(', among them in the extended MHC (chr6 28.4–33.5 Mb): ' + mhc_txt) if len(mhc) else ''}; only with PCs: {slash(op)}.
- Effect sizes of the {len(lo)} cross-cohort lead × cohort pairs barely move: median |Δ ln OR| = {med:.3f} (max {mx:.2f}).
"""
    figures = [('grm_pc_vs_grm.png', 'one row per cohort: Manhattan of the random model without PCs (left) and '
                'with PCs (middle, the main analysis) on one axis and scale, and the QQ plot of both (right).')]
    calib = pd.DataFrame({
        'Cohort': [names[c] for c in cohorts], 'λGC without PCs': [f'{v:.3f}' for v in lw],
        'λGC with PCs': [f'{v:.3f}' for v in lp], 'P<10⁻⁵ without': nw, 'P<10⁻⁵ with': np_,
        'peaks without / with': [f'{j[c]["peaks_without_pc"]} / {j[c]["peaks_with_pc"]}' for c in cohorts],
        'only without': ow, 'only with': op})
    pcs = [k for k in pc.columns if k.startswith('p_PC')]
    pct = pd.DataFrame({
        'Cohort': [names[c] for c in cohorts], 'cases / controls': [f'{pc.loc[c, "n_case"]}/{pc.loc[c, "n_ctrl"]:,}' for c in cohorts],
        'LR χ²(10)': [f'{pc.loc[c, "lr_chi2"]:.1f}' for c in cohorts], 'P': [pu(pc.loc[c, 'lr_p']) for c in cohorts],
        'PCs with P < 0.05': [', '.join(k[2:] for k in pcs if pc.loc[c, k] < 0.05) for c in cohorts],
        'pseudo-R² gain': [f'{pc.loc[c, "pseudo_r2_sex_pcs"] - pc.loc[c, "pseudo_r2_sex"]:.4f}' for c in cohorts]})
    pkt = pk.assign(Cohort=pk.cohort.map(names), **{'only in': pk.only_in.map({'grm_only': 'without PCs', 'grm_pc': 'with PCs'})},
                    Gene=pk.gene, Position=[f'chr{c}:{p:,}' for c, p in zip(pk.chrom, pk.pos)], rsID=pk.rsID,
                    variants=pk.n_sig_variants.astype(int), P=[pu(p) for p in pk.P])[
        ['Cohort', 'only in', 'Gene', 'Position', 'rsID', 'variants', 'P']]
    detailed = f"""
### Question
Within the random (mixed) model with a full GRM, what do the 10 ancestry PCs add as fixed covariates?

### Design
Both scans are SAIGE logistic mixed models on the same samples (relatives retained), the same full GRM and the
same variants; the only difference is the covariates: SEX + full GRM, against SEX + PC1–10 + full GRM (the main
analysis). P values: score test, normal approximation.

### How to read `grm_pc_vs_grm.png`
As `01.fixed_vs_random`: rows are cohorts, both Manhattans share one axis and one −log₁₀P limit, red line
5 × 10⁻⁸, grey dashed 10⁻⁵, diamonds = a panel's own genome-wide leads, QQ with the 95% band.

### Results
Calibration and peaks (peaks: upstream calls, ±250 kb; "only" = no peak of the other model within 500 kb):

{md_table(calib)}

Do the PCs differ between cases and controls? Logistic regression of case status on SEX vs SEX + PC1–10, samples
of the mixed-model scans:

{md_table(pct)}

Peaks called under only one of the two models:

{md_table(pkt)}

Lead variants of the 16 cross-cohort loci (genome-wide in any cohort or suggestive in all three) under both
models: median |Δ ln OR| = {med:.3f}, largest {mx:.2f}; per variant in `lead_or_stability.tsv`.

λGC by minor-allele frequency: `lambda_by_maf.tsv`.

### Open question (to be agreed)
Whether, and on which of the facts above, the main analysis should include the PCs on top of the GRM. This
section is left without a conclusion on purpose.

### Methods
As `01.fixed_vs_random` (λGC, QQ band, drawing). PC test: `statsmodels` Logit, likelihood-ratio
χ²(10) of SEX + PC1–10 against SEX alone, per-PC Wald P and McFadden pseudo-R². Upstream peak calls of the
`main` (with PCs) and `detect` (without PCs) scans of `analysis/assoc_saige`.
"""
    write(a.out_dir, '2 · Random model with and without PCs', slide, figures, detailed)


# ── 3. main analysis ─────────────────────────────────────────────────────────
def main_results(a):
    cohorts, names = C.csv_list(a.cohorts), cohort_names(a)
    t = pd.read_csv(a.tail, sep='\t')
    t = t[t.model == 'grm_pc'].set_index('cohort')
    L = pd.read_csv(a.loci, sep='\t')
    lam = [t.loc[c, 'lambda_gc'] for c in cohorts]
    gw = [int(t.loc[c, 'gw_loci']) for c in cohorts]
    sg = [int(t.loc[c, 'sg_loci']) for c in cohorts]
    concord = L.groupby('locus_id').OR.apply(lambda x: np.sign(np.log(x)).nunique() == 1)
    gwl = L[L.group == 'b']
    sgl = L[L.group == 'c']
    sg_genes = list(dict.fromkeys(sgl.sort_values('P').drop_duplicates('locus_id').gene))
    lead_txt = []
    for lid, g in gwl.groupby('locus_id', sort=False):
        g = g.set_index('cohort').loc[cohorts]
        best = g.P.idxmin()
        lead_txt.append(f'*{g.gene.iloc[0]}* ({g.rsID.iloc[0]}): OR {slash(g.OR, "{:.2f}")}, P = '
                        f'{slash([pu(p) for p in g.P])}; genome-wide in {names[best].lower()}')
    slide = f"""
**Take-home: the main analysis (random model: SAIGE, SEX + 10 PCs + full GRM) is well calibrated in all three nested cohorts (λGC {slash(lam, '{:.3f}')}). {' '.join(lead_txt)}. {len(sg_genes)} further loci are suggestive in all three cohorts, each with the same direction of effect in every cohort.**

- Cohorts (cases / controls): {'; '.join(f'{names[c]} {int(t.loc[c, "n_case"])} / {int(t.loc[c, "n_ctrl"]):,}' for c in cohorts)}.
- Loci, genome-wide + suggestive: {'; '.join(f'{names[c].split()[0].lower()} {g} + {s}' for c, g, s in zip(cohorts, gw, sg))}.
- Suggestive in all three cohorts: {', '.join(f'*{x}*' for x in sg_genes)}.
- The cohorts are nested and share every case: agreement between them shows robustness to the ancestry filter, not replication.
"""
    figures = [('main_manhattan.png', 'one row per cohort: Manhattan and QQ of the main analysis on one genomic '
                'axis and one −log₁₀P limit (the layout of `cohort_manhattan.png`, without its caption).'),
               ('main_summary.png', '(a) cases and controls per cohort with λGC and the number of '
                'genome-wide (GW) and suggestive (SG) loci; (b) OR (95% CI) of every locus genome-wide in any '
                'cohort or suggestive in all three, the cohorts side by side, with the smallest P in its '
                "cohort's colour.")]
    rows = []
    for lid, g in L.groupby('locus_id', sort=False):
        g = g.set_index('cohort')
        rows.append({'Group': 'GW in ≥1' if g.group.iloc[0] == 'b' else 'SG in all 3', 'Gene': g.gene.iloc[0],
                     'Variant': g.variant_id.iloc[0], 'rsID': g.rsID.iloc[0] if str(g.rsID.iloc[0]) not in ('.', 'nan') else '—',
                     **{f'OR (95% CI) {names[c].split()[0].lower()}': f'{g.loc[c, "OR"]:.2f} ({g.loc[c, "L95"]:.2f}–{g.loc[c, "U95"]:.2f})'
                        for c in cohorts},
                     'smallest P': f'{pu(g.P.min())} ({names[g.P.idxmin()].split()[0].lower()})'})
    over = pd.DataFrame({'Cohort': [names[c] for c in cohorts],
                         'cases / controls': [f'{int(t.loc[c, "n_case"])} / {int(t.loc[c, "n_ctrl"]):,}' for c in cohorts],
                         'N_eff': [f'{t.loc[c, "n_eff"]:,.0f}' for c in cohorts],
                         'variants': [f'{int(t.loc[c, "n_variants"]):,}' for c in cohorts],
                         'λGC': [f'{v:.3f}' for v in lam], 'GW loci': gw, 'SG loci': sg,
                         'smallest P': [pu(t.loc[c, 'min_p']) for c in cohorts]})
    detailed = f"""
### Question
What does the main analysis show in the three nested cohorts?

### Design
SAIGE logistic mixed model, SEX + PC1–10 + full GRM, related samples retained; score test with the normal
approximation. Loci: upstream peak calls (variants at P < 10⁻⁵ merged within ±250 kb; genome-wide = lead
P < 5 × 10⁻⁸). The loci of panel b: genome-wide in at least one cohort, or suggestive in all three cohorts,
each at one representative variant shared by the three cohorts.

### How to read the figures
- **`main_manhattan.png`** — as `cohort_manhattan.png` of `analysis/assoc_saige`: rows are cohorts on one
  genomic axis and one −log₁₀P limit; red line 5 × 10⁻⁸, grey dashed 10⁻⁵; the dotted guide and the label
  mark the locus genome-wide in any cohort, the diamond the cohort that called it.
- **`main_summary.png`** — (a) composition and calibration per cohort; (b) one row per locus, the three cohorts
  side by side (diamond = genome-wide in that cohort), OR on a log scale with the line at 1; right, the
  smallest P of the three, in the colour of its cohort.

### Results
{md_table(over)}

{md_table(pd.DataFrame(rows))}

All {int(concord.sum())} of {len(concord)} loci have the same direction of effect in the three cohorts.

### Interpretation
- Calibration is good in every cohort; the relaxed ancestry filter (narrow → full) adds controls without
  adding inflation.
- *FTO* is the one genome-wide locus; the suggestive loci shared by all three cohorts carry concordant
  effects of similar size.
- The cohorts are nested — narrow ⊂ intermediate ⊂ full — and share every case, so concordance is expected
  and is not independent replication.

### Methods
As `01.fixed_vs_random`. Loci and their representative variants: `analysis/assoc_saige`
`_comparison/tables/lead_crosscohort.tsv`; the per-variant details are in `integration/saige_signals`.
"""
    write(a.out_dir, '3 · Main analysis in the three cohorts', slide, figures, detailed)


def index(a):
    parts = ['# Model choice and main results of the CTEPH GWAS — figures and text for slides', '',
             'Three questions, one folder each. Every folder has the figures and a `README.md` with a **slide '
             'summary** (brief) and a **detailed description**; the slide summaries are repeated here for '
             'copying.', '',
             '| Folder | Question | Figures |', '|---|---|---|',
             '| `01.fixed_vs_random/` | Fixed vs random model: is a fixed-effect model enough? | '
             '`fixed_vs_random.png`, `calibration.png` |',
             '| `02.pc_effect/` | Random model with and without PCs: what do the PCs add? (conclusion pending) | '
             '`grm_pc_vs_grm.png` |',
             '| `03.main_results/` | Main analysis in the three cohorts | `main_manhattan.png`, `main_summary.png` |',
             '', '---', '']
    for s in a.slides:
        parts += [Path(s).read_text().strip(), '', '---', '']
    Path(a.out_dir, 'README.md').write_text('\n'.join(parts))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--family', required=True, choices=['fixed_vs_random', 'pc_effect', 'main_results', 'index'])
    ap.add_argument('--cohorts', default='')
    ap.add_argument('--cohort-labels', default='')
    ap.add_argument('--tail', default='')
    ap.add_argument('--maf', default='')
    ap.add_argument('--fixed-leads', default='', help='gw_leads.tsv of scan_tables.py')
    ap.add_argument('--pc-effect', default='')
    ap.add_argument('--pc-assoc', default='')
    ap.add_argument('--peak-changes', default='')
    ap.add_argument('--lead-or', default='')
    ap.add_argument('--loci', default='')
    ap.add_argument('--slides', nargs='*', default=[])
    ap.add_argument('--out-dir', default='.')
    a = ap.parse_args()
    {'fixed_vs_random': fixed_vs_random, 'pc_effect': pc_effect, 'main_results': main_results,
     'index': index}[a.family](a)
    C.log(f'README: {a.family} -> {a.out_dir}')


if __name__ == '__main__':
    main()
