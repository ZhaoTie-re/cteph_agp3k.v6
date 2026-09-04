#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : ONE document per figure family.
#
#           The variant-set figure is drawn once per cohort and the gene-scan
#           figure once per cohort x stratum. Their prose -- what each panel
#           shows, how to read it, what it cannot establish -- is a property of
#           the figure, identical for every member; only the numbers differ. A
#           sidecar per PNG would be the same page copied six times. So each
#           plotting task drops its numbers beside its PNG as <name>.stats.json
#           (figure_doc.write_stats) and this step writes the family's README.md
#           once, with every member as one row of one table -- and the members
#           become comparable down a column, which separate files never are.
#
#           The prose lives here and not in the plotting scripts because it
#           describes the family, not any one rendering.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_shared' / 'scripts'))
import figure_doc                      # noqa: E402

METHODS_REF = '../../../docs/METHODS.md'     # results/figures/<family>/README.md -> docs/


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family', required=True, choices=sorted(FAMILIES))
    p.add_argument('--stats', nargs='+', required=True,
                   help='the members\' <name>.stats.json files')
    p.add_argument('--figure-dir', required=True,
                   help='the family directory name, e.g. 02.gene_scan (for the subject line)')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--min-num-var', type=int, default=2)
    p.add_argument('--out', default='README.md')
    return p.parse_args()


def read_members(paths):
    """One row per member; columns in first-seen order of the value labels."""
    rows, order = [], []
    for p in paths:
        d = json.load(open(p))
        row = {'figure': d['peak']}
        for label, value in d['values']:
            if label not in order:
                order.append(label)
            row[label] = value
        rows.append(row)
    df = pd.DataFrame(rows, columns=['figure'] + order)
    return df.sort_values('figure').reset_index(drop=True)


# ── family prose ─────────────────────────────────────────────────────────────
def _variant_sets(a, df):
    return dict(
        title='Variant sets — what the map hands to the tests',
        question=('Of the variants that pass the cohort\'s minor-allele-count floor, which enter '
                  'a gene set, under which consequence, and how large are the sets each stratum '
                  'produces?'),
        panels=[
            ('a', 'Call set by snpEff impact',
             'Every minAC-passing variant of the cohort, once, at its most severe annotation. '
             'MODIFIER (intergenic, deep intronic) is the vast majority and never enters a gene '
             'set; the three classes that can are HIGH, MODERATE and LOW.'),
            ('b', 'Consequences by impact class',
             'The most frequent snpEff consequence terms within HIGH, MODERATE and LOW, each '
             'term counted once per variant at its most severe annotation; the remaining terms '
             'of each class are folded into an "Other" row whose count of distinct terms is '
             'given. Compound terms are shown by their first component with the number of '
             'additional components in brackets; the full strings are in '
             'annotation_summary.tsv.'),
            ('c', 'Gene sets and variants per stratum',
             'For each impact stratum, the number of gene sets the map assigns (any size) '
             'beside the number that are tested (at least MinNumVar mapped variants), and the '
             'number of distinct tested variants. The tested count is the Bonferroni '
             'denominator of that cohort x stratum.'),
            ('d', 'Set-size distribution',
             'The fraction of tested sets with at least x variants, per stratum, on a log axis. '
             'Most tested sets are small; the dashed guide is MinNumVar, below which a set is '
             'not tested and not counted.'),
        ],
        interpretation=(
            'The figure shows the experiment that was designed, before any association was '
            'run: which variants qualify, how they distribute over consequence classes, and how '
            'many genes each stratum can test. Because MinNumVar is applied before the '
            'denominator is fixed, a set absent here is absent from the family the thresholds '
            'correct over. The three cohorts differ almost entirely in controls, so their maps '
            'differ mainly by which rare variants clear the minor-allele-count floor.'),
        numbers=[('members', len(df)), ('alpha', a.alpha), ('MinNumVar', a.min_num_var)],
        tables=[('Per-cohort values', df)],
        reading=[
            'Read (c) first: the tested-set count per stratum is the denominator every '
            'p-value in that cohort x stratum is corrected over.',
            'Then (a) and (b): which consequences actually make up the sets — a stratum '
            'dominated by one term is a different experiment from a mixed one.',
            '(d) says how small the sets are; a gene with two variants is tested, but its '
            'statistic is a count over a handful of alleles (METHODS §4.4).',
        ],
        limits=[
            'It shows annotation, not association: nothing here says which gene carries '
            'signal.',
            'Counts are per cohort and depend on that cohort\'s minor-allele-count floor; they '
            'are not comparable across cohorts as if the callset were fixed.',
        ],
        defs=[],
        model='snpEff canonical map, per-gene most-severe consequence; sets need >= MinNumVar variants',
        methods_ref=METHODS_REF)


def _gene_scan(a, df):
    return dict(
        title='Gene-based scans — CMC burden and SKAT-O, every cohort x stratum',
        question=('Where does each scan show association, are the two statistics calibrated, '
                  'and do they agree about the genes they single out?'),
        panels=[
            ('a', 'CMC burden across the genome',
             'One point per gene at the midpoint of its mapped variant span, chromosomes '
             'banded, y = -log10 of the CMC score-test P (logistic model, sex + ancestry PCs, '
             'minor alleles collapsed to carrier / non-carrier). The solid line is the '
             'Bonferroni threshold alpha / n_genes_mapped fixed at the map barrier; the dashed '
             'line is the Benjamini-Hochberg cut of THIS statistic — the largest P '
             'it rejected — over the same gene family; where BH rejected nothing it has no cut, '
             'and the panel says so instead of leaving the reader to notice a missing line. '
             'Genes called by both rules are drawn in '
             'the full accent, genes called by BH only in the lighter one, and EVERY called '
             'gene -- and only a called gene -- is named above the panel in the same two '
             'tones, written vertically so that every panel of the family is labelled the '
             'same way. A scan that called nothing names nothing.'),
            ('b', 'SKAT-O across the genome',
             'The same construction over the SKAT-O P-values (same model, same sets, same '
             'samples, mean-imputed dosages for the kernel). The panel shares (a)\'s x mapping, '
             'y-limit and Bonferroni line, so a taller point is a stronger point; the BH cut is '
             'this statistic\'s own. Both Manhattan data boxes are the same fixed height in every '
             'member of the family, and the y limit is the largest -log10 P this variant set '
             'reached in ANY cohort (scan_scale.tsv), so a point may be compared across the three '
             'cohorts of one set; the threshold lines stay each cohort\'s own.'),
            ('c', 'Calibration',
             'Observed against expected -log10 P for both statistics over the genes each '
             'returned, with the pointwise Beta(i, n-i+1) 95 % band of the null. The legend '
             'carries the genomic-control inflation factor lambda_GC of each statistic, read '
             'from scan_qc.tsv. The split of lambda by set size (n_var = 2 against >= 3) is in '
             'the calibration table below, not in the figure: for two-variant genes the '
             'normal approximation of a discrete statistic is anti-conservative near P = 0.5, '
             'which inflates a median-based lambda without touching the tail where calls are '
             'made (METHODS §4.4).'),
            ('d', 'CMC against SKAT-O, gene by gene',
             'Genes returned by both statistics, matched on set name. A gene is accented when '
             'BH calls it under either statistic; the solid lines are the Bonferroni threshold '
             'on both axes, the dashed lines each statistic\'s BH cut. A gene far above the '
             'diagonal carries a signal the kernel sees and the collapse does not (a subset of '
             'variants, or opposite directions); a gene below it is one where every minor '
             'allele points the same way. No gene is named here -- the names are in (a) and (b), '
             'against the genome and the thresholds; this panel marks the called genes in red '
             'and counts them.'),
        ],
        interpretation=(
            'Each cohort x stratum is one scan judged by one threshold, alpha / n_genes_mapped, '
            'fixed before any P existed, and by BH q < alpha over the same family. A gene '
            'called by one statistic and not the other is therefore a statement about the '
            'statistics, not about two thresholds. Calibration bounds everything: a curve '
            'leaving the band early is producing small P-values everywhere. None of this makes '
            'a gene a CTEPH gene — cases and controls were sequenced on different platforms '
            'with no overlap, so every P carries that confound and a called gene is a candidate. '
            'Which candidates survive the change of statistic and the addition of controls is '
            'the question of the robust-genes figure.'),
        numbers=[('members (cohort x stratum)', len(df)), ('alpha', a.alpha),
                 ('MinNumVar', a.min_num_var)],
        tables=[('Per-scan values', df)],
        reading=[
            'Read (a) and (b) against the same solid line, then against each other.',
            'Go to (c) before believing either: lambda_GC is reported, not targeted, and the '
            'table gives its set-size split.',
            'Use (d) to see which genes the two statistics disagree about; those are the '
            'genes whose variant-level structure matters.',
            'A gene-set P says nothing about which variant carries it; the per-gene figures do.',
        ],
        limits=[
            'It does not establish that any gene is associated with CTEPH (platform confound).',
            'Point HEIGHTS are comparable across the three cohorts of one variant set, whose '
            'y axis is shared (scan_scale.tsv), and across (a) and (b) within a figure. They '
            'are NOT comparable between variant sets, and no DECISION is comparable anywhere: '
            'each scan has its own denominator, so the threshold lines sit at different '
            'heights in different scans.',
            'The diagonal in (d) is not a null expectation; two statistics on the same data are '
            'correlated by construction.',
        ],
        defs=[],
        model=('rvtest --burden cmc,cmcWald and --kernel skato, logistic, sex + ancestry PCs, '
               '--impute mean; threshold = alpha / n_genes_mapped per cohort x stratum; '
               'called = BH q < alpha over the same family'),
        methods_ref=METHODS_REF)


FAMILIES = {'variant_sets': _variant_sets, 'gene_scan': _gene_scan}


def main():
    a = parse_args()
    df = read_members(a.stats)
    kw = FAMILIES[a.family](a, df)
    subject = (f'**Figures:** `{a.figure_dir}/` — {len(df)} figure(s): '
               + ', '.join(f'`{f}`' for f in df['figure']) + '. One document covers the family; '
               'the numbers of each member are one row of the table below.')
    figure_doc.write_doc(None, out_path=a.out, subject=subject, **kw)
    print(f'[catalogue] {a.family}: {len(df)} member(s) -> {a.out}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:                       # noqa: BLE001
        print(f'Error in catalogue: {e}', file=sys.stderr)
        sys.exit(1)
