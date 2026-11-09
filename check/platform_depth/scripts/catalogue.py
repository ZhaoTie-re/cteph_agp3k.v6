#!/usr/bin/env python3
"""One README.md per figure family: what each figure shows and where its numbers come from.

The descriptions are fixed here; the numbers quoted come from each figure's <name>.stats.json, so a
README never states a value the figure was not drawn from.
"""
import argparse
import json
from pathlib import Path

FAMILIES = {
    '01.depth': ('Depth and the down-sampling target', {
        'depth_platforms': 'a, effective depth (MQ ≥ 20, BQ ≥ 10, duplicates excluded) of every full-mainland genome by '
                           'platform, median as a bar; the dashed line is the target D*, the median of the HiSeq X '
                           'controls; arrows take every down-sampled platform\'s median to it (all case platforms but '
                           'DNBseq-G400RS 15x). b, the fraction of read pairs each down-sampled case keeps: the deep '
                           'cases each to D*, the DNBSeq-T7 cases one platform fraction. c, the read-length audit: '
                           'measured read length against the length the info sheet\'s Observed_Depth implies '
                           '(DNBSeq-T7 is a 100 bp library).',
        'depth_local': 'Median effective depth at every calling interval: HiSeq X controls, DNBSeq-T7 cases before and '
                       'after their platform fraction, and the deep cases after down-sampling (their depth there × '
                       'their fraction).',
    }),
    '02.downsampling': ('Down-sampling: genotypes and signals', {
        'concordance': 'Genotype concordance among genotypes called in both, heterozygotes retained and call rate of '
                       'the down-sampled cases: re-processing (baseline vs original) and the depth effect (each '
                       'replicate vs baseline), by variant class.',
        'concordance_platforms': 'The same metrics over every fixed site, per platform of the down-sampled cases '
                                 '(DNBSeq-T7, NovaSeq, DNBSeq-G400RS 30x).',
        'gwas_leads': '-log10 P of each locus lead per cohort: original genotypes (open), the f = 1 baseline (square) and '
                      'the range over replicates (line); dotted 10^-5, dashed 5 × 10^-8.',
        'hla_dqb1_57v': 'DQB1:57:V, P divided by the scan\'s GEC threshold (below 1 = significant): original, baseline, '
                        'replicate range, per cohort.',
    }),
    '03.platform': ('Platform sensitivity', {
        'forest.narrow_mainland': 'Odds ratio of each GWAS lead and DQB1:57:V in each case-platform stratum (that '
                                  'platform\'s cases vs every control), the pooled estimate (diamond), and Q and I² '
                                  'across the strata; narrow mainland.',
        'forest.intermediate_mainland': 'As above, intermediate mainland.',
        'forest.full_mainland': 'As above, full mainland.',
    }),
    '04.array': ('SNP-array evidence, per cohort', {
        'array_evidence.narrow_mainland': 'The unrelated narrow-mainland genomes with SNP arrays (DNBSeq-T7 cases on '
                                          'ASA24v1, HiSeq X controls on JSA24v1), WGS against the array data in the '
                                          'same people. a, OR on WGS against OR on the array data at one marker per '
                                          'signal: the lead typed, else imputed (R² ≥ 0.8; DQB1:57:V imputed with '
                                          'HIBAG), else an LD proxy (r² ≥ 0.8); cases and controls are read from the '
                                          'same data; the line is identity.',
        'array_evidence.intermediate_mainland': 'As above, intermediate mainland.',
        'array_evidence.full_mainland': 'As above, full mainland (the people of check/array_crosscheck; gate G10), '
                                        'with b, the driver variants of the gene sets (full mainland, which contains the '
                                        'other two cohorts): WGS case carriers and those the imputed arrays confirm '
                                        '(dosage ≥ 0.5); controls are in the imputation panel and not shown.',
    }),
    '05.controls': ('Control sources', {
        'sources.narrow_mainland': 'Odds ratio of each GWAS lead and DQB1:57:V with every case against the controls of '
                                   'one source (Nagahama, BBJ, HTLV-1, ACC, lung cancer; each with its own SAIGE null '
                                   'model), all controls (diamond), and Q and I² across the sources (they share the '
                                   'cases); narrow mainland.',
        'sources.intermediate_mainland': 'As above, intermediate mainland.',
        'sources.full_mainland': 'As above, full mainland.',
    }),
    '06.evidence': ('The evidence matrix and the follow-up of its flags', {
        'evidence_matrix': 'Every signal (rows: GWAS leads, DQB1:57:V, position DQB1 57, Tier 1–2 genes in their M/H '
                           'sets, the main analysis) against every check (columns), read by the rules of 14.evidence/evidence.rules.tsv: '
                           'blue consistent, red flag, grey not testable or not checkable, white not '
                           'applicable. Three cells per check = narrow / intermediate / full mainland. A check flags '
                           'a signal when its rule holds in any one cohort; nothing is removed.',
        'z_observed_expected': 'Every leave-one-out estimate: z observed (sign set so that the pooled direction is '
                               'positive) against |z| expected = pooled β / SE without the left-out group. a–c, leave '
                               'one case platform out; d–f, leave one control source out; per cohort. The line is '
                               'identity; the shaded corner (|z| expected ≥ 2.8, z observed < 1.96) is the flag rule. '
                               'Genes use the CMC burden β.',
        'gene_carriers.narrow_mainland': 'Tier 1–2 genes in their tiered sets, narrow mainland: carrier share per '
                                         'case platform and per control source (counts beside each point), the pooled '
                                         'groups, and the share ToMMo 61KJPN predicts from its allele frequencies '
                                         '(diamond). "set missing": genotypes the cohort filters set missing although '
                                         'the pre-QC record carries the ALT.',
        'gene_carriers.intermediate_mainland': 'As above, intermediate mainland.',
        'gene_carriers.full_mainland': 'As above, full mainland.',
        'lead_frequencies.narrow_mainland': 'GWAS leads, narrow mainland: ALT allele frequency per case platform and '
                                            'control source, the pooled groups (dotted lines: all cases red, all '
                                            'controls blue), ToMMo 61KJPN and 1KG JPT.',
        'lead_frequencies.intermediate_mainland': 'As above, intermediate mainland.',
        'lead_frequencies.full_mainland': 'As above, full mainland.',
        'v57_replicates': 'The deep cases whose DQB1:57:V dosage differs between any two arms (original typing, the '
                          'f = 1 baseline, the replicates): the DQB1 call, the dosage and the HLA-HD confidence; '
                          'red one V copy more than the baseline, blue one fewer. The deep cases belong to all three '
                          'cohorts and their typings are the same in each, so the figure is drawn once.',
    }),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--figures', required=True)
    a = ap.parse_args()
    root = Path(a.figures)
    for fam, (title, figs) in FAMILIES.items():
        d = root / fam
        if not d.exists():
            continue
        lines = [f'# {title}', '', f'Figures of `check/platform_depth`, family `{fam}`. Every number drawn is read from '
                 'the result tables (docs/OUTPUTS.md); `<name>.stats.json` holds the values each figure was drawn from.', '']
        for name, text in figs.items():
            png = d / f'{name}.png'
            if not png.exists():
                continue
            st = d / f'{name}.stats.json'
            lines += [f'## {name}', '', f'![{name}]({name}.png)', '', text, '']
            if st.exists():
                lines += ['Drawn from: `' + json.dumps(json.loads(st.read_text()), default=str)[:400] + '`', '']
        (d / 'README.md').write_text('\n'.join(lines))
        print(f'[catalogue] {fam}: README written')


if __name__ == '__main__':
    main()
