# Typing accuracy differs between the groups; it is not depth and not the panel

**Figure file:** `typing_confound.png`

## The question this figure answers

The controls are typed measurably worse than the cases. Is that sequencing depth — which could be controlled for — or an artefact of the reference panel, or something that can be neither corrected nor explained away?

## Panels

**(a) Unconfirmed share against measured depth**

One point per sequencing platform, x its median measured depth, y the share of called chromosomes carrying an allele the reference panel never sees; area is proportional to √n. Read it twice: along the case platforms alone, phenotype is held fixed and the share falls as depth rises — depth is real. Then compare the control platform with the case platform beside it at the same depth. The note in the panel gives the two group medians and a Mann–Whitney test of their equality.

**(b) The comparison at matched depth**

The two platforms carrying the comparison, restricted to samples between 17× and 21× so their median depths are equal. If depth were the cause the bars would be level. Fisher exact on the 2×2 of unconfirmed against confirmed chromosomes.

**(c) Per sample**

How many of a sample's called chromosomes at the reference loci sit on an unconfirmed allele. Shown as a distribution because a mean hides that most samples have none and a tail has several. The denominator is printed on the axis and is read from the data, not assumed.

**(d) Both reference panels**

The control and case shares computed twice from the same calls, once against each panel, with the control/case ratio printed above each pair. This is the panel that licenses the reading of every other one: if the ratio moved with the panel, the finding would be about the reference rather than about the typing.

## Interpretation

Three explanations are available and two of them are ruled out here. Depth genuinely affects typing — panel (a) shows it cleanly within the case platforms, where phenotype cannot be responsible — but cases and controls sit at the same measured depth and the gap survives matching on it, so that effect has nothing to act on between the groups. The reference panel genuinely sets the level — panel (d) shows the control level move from 4.63% against jMorp 61KJPN-HLA to 7.42% against 1000 Genomes JPT — but the ratio between the groups is stable across both (1.31 and 1.37), so the contrast is not an artefact of what the panel happens to carry. What is left is platform and cohort: every control is an AGP3K sample on HiSeqX, every case is a PH sample on something else. Those are perfectly confounded with each other and with phenotype, so which of read length, library chemistry, alignment reference or batch is responsible cannot be determined from these data. A common allele depleted more in controls than in cases reads as enrichment in cases — a false risk association from the typing alone.

## Values in this rendering

| quantity | value |
|---|---|
| controls, jMorp 61KJPN-HLA | 4.63% |
| cases, jMorp 61KJPN-HLA | 3.53% |
| ratio, primary panel | 1.31 |
| controls, 1000 Genomes JPT | 7.42% |
| cases, 1000 Genomes JPT | 5.40% |
| ratio, secondary panel | 1.37 |
| Fisher P, primary panel | 1.89e-07 |
| Fisher P, secondary panel | 7.61e-07 |
| control median depth | 18.72× |
| case median depth | 18.65× |
| Mann-Whitney P, depth | 0.194 |
| Fisher P, matched 17-21x | 3.58e-07 |

## Full statistics

**Every stratum behind the four panels, as written to typing_confound.tsv**

| stratum | n | median_depth | n_chr | n_unconfirmed | share_unconfirmed | fisher_p_matched | mannwhitney_p_depth | fisher_p_group |
|---|---|---|---|---|---|---|---|---|
| controls | 2,662 | 18.72 | 64,130 | 2,972 | 0.04634 |  |  |  |
| cases | 439 | 18.65 | 10,564 | 373 | 0.03531 |  | 0.194 | 1.89e-07 |
| platform:DNBseq-G400RS 15x | 31 | 14.17 | 732 | 78 | 0.1066 |  |  |  |
| platform:DNBSeq-T7 30x | 319 | 18.43 | 7,684 | 258 | 0.03358 |  |  |  |
| platform:HiSeqX 15x | 2,662 | 18.72 | 64,130 | 2,972 | 0.04634 |  |  |  |
| platform:NovaSeq 30x | 50 | 30.93 | 1,208 | 16 | 0.01325 |  |  |  |
| platform:DNBSeq-G400RS 30x | 39 | 36.25 | 940 | 21 | 0.02234 |  |  |  |
| matched 17-21x:HiSeqX 15x | 1,587 | 18.58 | 38,236 | 1,714 | 0.04483 |  |  |  |
| matched 17-21x:DNBSeq-T7 30x | 281 | 18.53 | 6,764 | 214 | 0.03164 | 3.58e-07 |  |  |
| 1000 Genomes JPT:controls | 2,662 | 18.72 | 26,618 | 1,974 | 0.07416 |  |  |  |
| 1000 Genomes JPT:cases | 439 | 18.65 | 4,388 | 237 | 0.05401 |  |  | 7.61e-07 |

## What this figure does *not* establish

- The absolute unconfirmed share is not an error rate. An allele a panel never carries may still be real and rare; panel (d) measures how much that matters by changing the panel. Only the difference between groups, measured on the same panel, is interpretable.
- Only the loci a panel carries are visible. Loci absent from a panel are not scored against it at all, which is why the two panels contribute different numbers of chromosomes per sample and why each locus list travels with its panel.
- The cohort is ancestry-restricted, so an unconfirmed allele can no longer be explained by the sample not being mainland Japanese. That removes an explanation; it does not remove the platform confound, which is what remains.
- The figure identifies the confound; it does not correct it. No within-data control exists, which is why OPEN_QUESTIONS.md §2 — typing samples of known HLA type — is now the only remaining option.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
