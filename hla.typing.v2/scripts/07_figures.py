#!/usr/bin/env python3
"""Step 07 — the component's figures (journal width, paper type scale).

results/figures/
  01.typing_overview.png       call state, source of the 2-field allele, call confidence, per gene
  02.read_evidence.png         unique-read balance; DRB3/4/5 two-allele calls by DRB1-implied copies;
                               DRB1-implied vs HLA-HD copy number
  03.panel_frequencies.png     our controls vs ToMMo 61KJPN, one panel per locus
  04.locus_summary.png         the three tier criteria and the tier, per locus
  05.panel_ceiling.png         discordance with 1KG JPT: our controls vs ToMMo
  06.confidence_vs_panel.png   alleles absent from ToMMo, by call confidence
  07.sample_qc.png             read pairs and low-confidence calls per sample
  08.shift_direction.png       alleles >= 1 point apart from ToMMo: where 1KG JPT lies
The catalogue that explains each figure is results/figures/README.md.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import figlib as F  # noqa: E402


def main():
    d = F.load()
    o = C.D_FIG
    F.fig_overview(d, o / '01.typing_overview.png')
    F.fig_read_evidence(d, o / '02.read_evidence.png')
    F.fig_panel_frequencies(d, o / '03.panel_frequencies.png')
    F.fig_locus_summary(d, o / '04.locus_summary.png')
    F.fig_panel_ceiling(d, o / '05.panel_ceiling.png')
    F.fig_confidence_vs_panel(d, o / '06.confidence_vs_panel.png')
    F.fig_sample_qc(d, o / '07.sample_qc.png')
    F.fig_shift_direction(d, o / '08.shift_direction.png')
    print(f'[figures] 8 figures -> {o}')


if __name__ == '__main__':
    main()
