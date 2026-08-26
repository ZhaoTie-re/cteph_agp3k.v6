# HLA typing quality, by sequencing platform

**Figure file:** `typing_qc.png`

## The question this figure answers

Does typing quality differ between the platforms, and therefore between cases and controls, in a way that could imitate a finding?

## Panels

**(a) Call rate by platform**

The fraction of (sample, gene) pairs that produced a genotype, per sample, grouped by platform. A platform whose distribution sits lower is producing fewer calls, and since platform is confounded with phenotype here that would be indistinguishable from a real difference.

**(b) Field resolution by platform**

The mean number of IMGT fields resolved per call. Shorter reads carry less information for separating similar alleles, so a read-length effect appears here before it appears in the call rate.

**(c) Call rate by gene**

Which loci are typed reliably. The pseudogenes and the DRB paralogues are expected to be low; the classical loci are not.

**(d) Polymorphic residue positions by gene**

How many IMGT positions actually vary in this cohort — the number of testable sites each gene contributes.

## Interpretation

Read this by platform. Cases were sequenced on DNBSeq and NovaSeq, controls entirely on HiSeqX, with no overlap, so any case/control difference in typing quality is confounded by construction and cannot be attributed. Platform is the axis on which a technical gradient is visible as itself. Flat metrics across platforms are evidence that read length and depth are not driving the typing; a gradient is the size of the problem, stated.

## What this figure does *not* establish

- This figure says nothing about typing ACCURACY. It measures whether a call was produced and how deeply it resolved, not whether it is right. Accuracy needs an external truth set.
- A flat call rate does not rule out a systematic bias toward particular alleles, which would need the same truth set to detect.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
