# assoc_hla — this study's numbers

Everything here is measured from this run's outputs. Where a number moved during development,
the reason is given, because a number that changed silently is worse than one that never
existed.

## The design as run

| | narrow_mainland | intermediate_mainland | full_mainland |
|---|---|---|---|
| samples | 2,195 | 2,508 | 3,101 |
| cases / controls | 419 / 1,776 | 429 / 2,079 | 439 / 2,662 |
| HLA-typed | **100 %** | **100 %** | **100 %** |

The cohorts are **nested** — narrow ⊂ intermediate ⊂ full — so agreement between them is not
replication. It is robustness to where the ancestry boundary is drawn, and nothing more.

That all three are fully typed is worth stating plainly: the HLA results sit on exactly the
sample sets `assoc_plink2` and `assoc_saige` scanned, so a comparison against the genome-wide
signal is a comparison of loci, not of samples.

## What is tested

16 genes: A, B, C, DMA, DMB, DOA, DOB, DPA1, DPB1, DQA1, DQB1, DRA, DRB1, E, F, G.

- **DRB3, DRB4, DRB5 excluded** from the primary analysis for the three reasons in METHODS §4.
- **64 pseudogene allele columns excluded** — H, W, Y, T, J, K, L, V, DRB2, DRB6–9 and DPA2 have
  no GRCh38 primary-assembly locus. None carries residue data.

Marker counts, full cohort, after the frequency floor (MAC ≥ 20, call rate ≥ 0.95):

| | built | tested |
|---|---|---|
| allele | 2,175 | ~170 |
| residue | 2,693 | ~705 |
| amino-acid positions with ≥ 2 tested residues | 1,115 | ~305 |

## Three defects found by building this, and fixed

**1. An undetermined residue was being read as an absent one.** Before the fix, the strongest
residue signal in the study was `AA_A_175_G`, *P* = 2.5 × 10⁻⁸. HLA-A position 175 is 99.67 %
glycine, and of the 21 samples not homozygous for it **five — 24 % — had no determined residue
at all** rather than carrying alanine or arginine. Since typing quality is differential between
cases and controls here, that marker was measuring the typing artefact.

The fix masks a sample at every residue of a position whose dosages do not sum to 2 — 47,439
(sample, position) pairs on the full cohort — after which `AA_A_175_G` falls to MAC 16 and
leaves the scan through the ordinary frequency filter. Tested residues drop from 993 to 705.
**That reduction is the cost of the fix, and it is the right price.**

**2. A residue signal that was one non-existent allele.** After fix 1, the strongest residue
signal became `AA_A_178_T` at *P* = 8 × 10⁻⁹. Twenty chromosomes carry its minor state, 15 in
cases and 5 in controls — and **15 of those 20 carried `A*02:783`, a single allele that jMorp's
61,424 Japanese have never observed**. The allele layer had a flag for exactly this
(`in_reference`); the residue layer had none, and inherited the artefact without inheriting the
warning.

`frac_unconfirmed_allele` now carries it down: the share of a residue marker's minor-state
chromosomes whose allele is absent from the reference. Median 0.02 over 705 tested residues;
`AA_A_178_T` scores 0.85. It is a warning rather than a filter — a real Japanese allele the
panel never sampled scores high too — but no residue result should be reported without it.

**3. plink2 read every case as a control.** `pheno_covar.tsv` carries 0/1, which is what SAIGE
wants; plink2's own convention is 1 = control, 2 = case and it reads 0 as *missing*. Without
`--1` it loaded "419 controls" and no cases — and did not fail, it produced a clean run of
nonsense. Both engines disagreeing about phenotype coding is exactly the kind of thing that does
not announce itself.

## Significance: why three thresholds are not padding

`effective_tests.py` measures how many *independent* tests the scan really performed. On the
full cohort:

| class | tested | *M*<sub>eff</sub> | effective threshold | plain Bonferroni |
|---|---|---|---|---|
| allele | 170 | 123 (72 %) | 4.07 × 10⁻⁴ | 2.94 × 10⁻⁴ |
| residue | 993 | **149 (15 %)** | 3.36 × 10⁻⁴ | 5.04 × 10⁻⁵ |

The residue row is the argument. Nine hundred and ninety-three residue markers carry about
**149** independent tests, because the residues at one position are near-collinear by
construction and whole haplotypes travel together across genes. Plain Bonferroni would have
been **6.7× too strict**; 5 × 10⁻⁸ is not conservative here so much as mis-calibrated, since it
answers a question about a genome-wide scan of a million independent variants.

## λ_GC is above 1 here, and that is expected

The fixed model gives λ_GC ≈ 1.6–1.9 across cohorts. **This is not the usual diagnostic.**
λ_GC assumes the median marker is null. Inside 3.4 Mb of the most extreme LD in the genome,
with a real signal present, the median marker is correlated with the causal one and λ_GC is
inflated by construction. It is reported per cohort × model × class so the *comparison* between
models is readable, not so it can be used as a calibration test.

## Results

### Two amino-acid positions, consistent across every cohort and both models

The omnibus is where the result is. **The same two positions clear the threshold in all three
cohorts**, under both the effective-tests and the plain Bonferroni correction:

| position | df | narrow | intermediate | full |
|---|---|---|---|---|
| **DQB1:57** | 3 | 2.09 × 10⁻⁵ | 1.02 × 10⁻⁴ | 5.76 × 10⁻⁵ |
| **B:114** | 1 | 4.42 × 10⁻⁴ | 5.43 × 10⁻⁴ | 1.52 × 10⁻⁴ |

Both sit **inside the peptide-binding domain**, which is where a functional HLA signal is
expected and is drawn as a shaded band on `figures/02.omnibus/`. DQB1 position 57 is the
classical autoimmunity position. Neither was chosen in advance; both fell out of a scan of ~300
positions.

The residues carrying them, on the full cohort:

| marker | fixed | random (full GRM) | MAC | `frac_unconfirmed_allele` |
|---|---|---|---|---|
| `AA_DQB1_57_V` | OR 0.63, *P* 1.0 × 10⁻⁴ | OR 0.63, *P* 1.8 × 10⁻⁴ | 907 | 0.014 |
| `AA_B_114_N` | OR 1.35, *P* 1.3 × 10⁻⁴ | OR 1.34, *P* 3.5 × 10⁻⁴ | 2,425 | 0.067 |
| `AA_B_114_D` | OR 0.74, *P* 1.4 × 10⁻⁴ | OR 0.75, *P* 3.7 × 10⁻⁴ | 2,409 | 0.067 |

Common, fully determined, and clean on the artefact flag — the opposite profile from the two
signals the QC removed.

### What the two models add

Nothing changes when the GRM is fitted. Pooled log-OR correlation between the models is
**0.979**; the effect estimates are identical to two decimals; *P* values are marginally larger
under the GRM, as expected when a random effect absorbs some signal. λ_GC falls slightly
(allele class, full cohort: 1.60 → 1.47), so the GRM is doing something — it is just not doing
anything that changes a conclusion. **That is the useful result from running both**: it rules
out "this depends on the model" as an explanation.

λ_GC sits at 1.5–2.4 throughout and is **not** evidence of inflation here. It assumes the median
marker is null; inside 3.4 Mb of the strongest LD in the genome with a real signal present, the
median marker is correlated with the causal one. It is reported so the two models can be
compared, never applied as a correction.

### Nothing reaches 5 × 10⁻⁸

At the HLA-appropriate thresholds these positions are significant. At the genome-wide
convention none of them is. Both numbers are in every `significance.*.tsv`, and which one a
reader should use is argued in METHODS §6 rather than decided for them.

### The published Japanese CTEPH associations do not replicate

Screened against `info/prior_reports.tsv`, in the DVT-negative stratum the reports specified:

| allele | reported OR | ours (full, DVT−) | *P* | power for reported OR |
|---|---|---|---|---|
| HLA-B\*52:01 | 2.47 | **1.13 [0.77–1.65]** | 0.53 | ~100 % |
| HLA-DPB1\*02:02 | 5.07 | **0.89 [0.47–1.69]** | 0.72 | ~100 % |

**Both confidence intervals exclude the published point estimate**, in all three cohorts, and
the DVT-positive stratum is equally null (OR 1.02 and 0.98). This is a non-replication with the
power to be one, not an underpowered miss: the smallest odds ratios these strata could detect at
80 % power are 1.59 and 1.92.

Our control frequencies match theirs — B\*52:01 carriers 20.6 % against their 24 %,
DPB1\*02:02 8.1 % against their 6 % — so the alleles are being typed, they are simply not
enriched in cases. Against Kominami et al.'s 99 DVT-negative cases and 380 controls typed by
microsatellite, this is 145 cases against 2,662 controls typed directly from WGS.

## The limitation that outranks all of these

`hla.typing` OPEN_QUESTIONS §3: controls are typed measurably worse than cases, for a technical
reason perfectly confounded with phenotype, and a common allele depleted more in controls reads
as enrichment in cases. **No analysis in this component can separate that from a real effect.**
The `in_reference` flag marks the alleles most exposed to it — on the first narrow-cohort run the
single strongest allele signal, `HLA_F*01`, was flagged, being an unresolved one-field call on 28
chromosomes that a 61,424-person Japanese panel has never observed — but flagging is not
removing. Treat every allele result here as provisional until §2's accuracy measurement exists.
