# HLA typing quality, by sequencing platform

**Figure file:** `typing_qc.png`

## The question this figure answers

Does typing quality differ between the platforms, and therefore between cases and controls, in a way that could imitate a finding?

## Panels

**(a) Field resolution by platform**

The share of calls resolving to each IMGT field depth, with the 2-field share printed. Stacked rather than averaged: the per-platform MEAN sits between 2.84 and 2.95 for all five, a range nobody can act on, while the composition separates them cleanly — 10.3 % of calls stop at 2 fields on the weakest platform against 3.7 % on the strongest. No call in this cohort reached 4 fields, so that series is absent rather than drawn as a zero-width segment. Shorter reads carry less information for separating similar alleles, so a read-length effect would appear here first.

**(b) Ambiguity by platform**

Samples grouped by how many of their 33 genes HLA-HD emitted several candidate pairs for. It is the closest thing to a per-call confidence the tool provides, and it was previously computed but never plotted.

**(c) Outcome by gene**

The four outcomes, as a composition. `called` and `hemizygous` are both successes and their ratio is biology — DMA is hemizygous in 73 % of samples and HLA-A in 14 %, because DMA is nearly monomorphic in this population. `not_typed` means the gene is absent from that haplotype, which is why DRB3/4/5 look broken on a call-rate axis when `failed` is 0 for every gene in the cohort.

**(d) Polymorphic residue positions by gene**

How many IMGT positions actually vary in this cohort — the number of testable sites each gene contributes.

## Interpretation

Read this by platform. Cases were sequenced on DNBSeq and NovaSeq, controls entirely on HiSeqX, with no overlap, so any case/control difference in typing quality is confounded by construction and cannot be attributed. Platform is the axis on which a technical gradient is visible as itself. The gradient here is real but small — the platform ordering in (a) and (b) is the same one, and the weakest platform is the one sequenced shallowest. Its size is the size of the problem, stated. Whether it reaches the CALLS is a different question and a different figure: see typing_confound.png.

## What this figure does *not* establish

- Completeness is saturated on this cohort and therefore a weak quality signal: call rate takes five distinct values across the whole cohort and 97 % of samples sit on two of them. That is why every panel here is a composition and no panel is a mean — a mean over a saturated metric hides the only variation there is.
- This figure says nothing about typing ACCURACY. It measures whether a call was produced and how deeply it resolved, not whether it is right. Accuracy needs an external truth set.
- A flat call rate does not rule out a systematic bias toward particular alleles, which would need the same truth set to detect.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
