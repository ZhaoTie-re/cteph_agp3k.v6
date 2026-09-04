# Typing accuracy differs between the groups, and it is not depth

**Figure file:** `typing_confound.png`

## The question this figure answers

The controls are typed measurably worse than the cases. Is that sequencing depth — which could be controlled for — or something that cannot be?

## Panels

**(a) Unconfirmed share against measured depth**

One point per sequencing platform, x its median measured depth, y the share of called chromosomes carrying an allele the reference panel never sees; area is proportional to √n. Read it twice: along the case platforms alone, phenotype is held fixed and the share falls as depth rises — depth is real. Then compare the control platform with the case platform beside it at the same depth.

**(b) The comparison at matched depth**

The two platforms carrying the comparison, restricted to samples between 17× and 21× so their median depths are equal. If depth were the cause the bars would be level.

**(c) Per sample**

How many of a sample's ten called chromosomes at the five reference loci sit on an unconfirmed allele. Shown as a distribution because a mean hides that most samples have none and a tail has four or more.

**(d) Depth by group**

The measured-depth distributions of cases and controls with a Mann–Whitney test. This panel is why the downsampling control was retired: there is no depth difference to remove.

## Interpretation

Two effects are present and only one of them is a confounder. Depth genuinely affects typing — panel (a) shows it cleanly within the case platforms, where phenotype cannot be responsible. But cases and controls are at the same measured depth (panel d), so that effect has nothing to act on between the groups, and the gap survives matching on depth (panel b). What is left is platform and cohort: every control is an AGP3K sample on HiSeqX, every case is a PH sample on something else. Those are perfectly confounded with each other and with phenotype, so which of read length, library chemistry, alignment reference or batch is responsible cannot be determined from these data. A common allele depleted more in controls than in cases reads as enrichment in cases — a false risk association from the typing alone.

## What this figure does *not* establish

- The reference panel is 105 samples. An allele it never carries may still be real and rare, so the absolute level of the unconfirmed share is not an error rate. Only the difference between groups, measured on the same panel, is interpretable.
- Five loci only — A, B, C, DQB1, DRB1. DPA1, DPB1, DQA1 and DRB3/4/5 have no external reference and are invisible here.
- The figure identifies the confound; it does not correct it. No within-data control exists, which is why OPEN_QUESTIONS.md §2 — typing samples of known HLA type — is now the only remaining option.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
