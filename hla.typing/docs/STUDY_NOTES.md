# STUDY NOTES — hla.typing

This study's configuration and its measured numbers. Method and rationale are in
[METHODS.md](METHODS.md). Every number below is in `results/_run_info/facts.json` or in
[NUMBERS_ALLOWED.md](NUMBERS_ALLOWED.md); `verify.sh` §9 fails if one is not.

## Configuration

| | |
|---|---|
| **analysis cohort** | **3,101** — 2,662 AGP3K, 439 PH |
| how it was reached | `Flag_JHRPv6 = True` 3,592 → sample QC 3,569 (−23) → `full_mainland` 3,101 (−468) |
| excluded by sample QC | 23 heterozygosity outliers — HiSeqX 15x 18, DNBSeq-T7 30x 4, NovaSeq 30x 1 |
| excluded by ancestry | 468 — 455 AGP3K (all HiSeqX 15x), 13 PH; kept under `results/_superseded.3569/` |
| platforms in the cohort | HiSeqX 15x 2,662 · DNBSeq-T7 30x 319 · NovaSeq 30x 50 · DNBSeq-G400RS 30x 39 · DNBseq-G400RS 15x 31 |
| reference | `nagasaki_pipeline/data/hs38DH.fa` — 3,366 sequences, 525 `HLA*` contigs |
| extraction regions | 543 = MHC window + `chr6_*_alt` + `HLA*` + unmapped |
| typing | HLA-HD 1.7.1, `-m 100 -c 0.95` |
| dictionary | pinned IPD-IMGT/HLA 3.64.0, 46,005 alleles, `_N150` |
| residues | pinned IMGT protein alignments, same release, 19 genes |
| gene list | pinned `HLA_gene.split.3.50.0.txt`, 33 genes |
| primary frequency reference | jMorp 61KJPN-HLA — 61,424 individuals, 122,848 chromosomes per locus, 13 loci |
| secondary frequency reference | 1000 Genomes JPT — 105 individuals, 210 chromosomes, 5 loci |

**Both panels run on every execution**, each with its own locus list. See METHODS §10 for
why that is not redundancy, and why the two lists must never be shared.

## Sanity anchors — a refactor must not move these

The residue numbering is IMGT's, and these are how that is known:

| | |
|---|---|
| HLA-A positions 1–8 | **`GSHSMRYF`** — the canonical mature N-terminus |
| *DRB1* `*04:05:01` at 11/13/71/74 | **V / H / R / A** |
| *DRB1* `*15:01:01` at 11/13/71/74 | **P / R / A / A** |
| *DRB1* `*03:01:01` at 11/13/71/74 | **S / S / K / R** |
| `NXNE` in `residue_reference/` | **0 rows** — the signature of a fabricated residue |
| the frequency check against itself | feeding a panel in as if it were the calls gives `r = ρ = 1.000` and `max\|diff\| = 0.000` at every locus it carries |

These are properties of IMGT and of the code, not of the cohort, so they do not move when
the sample set does. That is what makes them anchors.

## Completeness, measured on the cohort

| | |
|---|---|
| (sample, gene) failures | **0** |
| gene calls | 48,684 called · 38,707 hemizygous · 14,942 not_typed · 0 failed |
| ambiguous (sample, gene) records | 2,558 |
| genes carrying residues | 19 of 33 |
| polymorphic residue positions | 1,297 of 5,957 reference positions |
| `allele_dosage.tsv` | 3,101 × 2,298 |
| `residue_dosage.tsv` | 3,101 × 3,125 |
| `residue_diplotype.tsv` | 3,101 × 1,297 |

Call rate is *k*/19 and is nearly saturated, which is why it is a weak quality metric here
and why every panel of `figures/typing_qc.png` is a composition rather than a mean.

## The confounding this component sits on

Platform is **perfectly separated from phenotype**: every control is HiSeqX, every case is
DNBSeq or NovaSeq. Nothing in this component can fix that, and no covariate can.

The two axes that could plausibly bias typing were checked before the component was built.

**Read length — measured, and the answer is no effect.** DNBSeq-T7 is a 100 bp library
against 150 bp everywhere else:

| platform | read length | *n* | call rate | field depth | median depth |
|---|---|---|---|---|---|
| DNBSeq-G400RS 30x | 150 bp | 39 | 0.919 | 2.9251 | 36.25× |
| **DNBSeq-T7 30x** | **100 bp** | **319** | **0.918** | **2.9307** | **18.43×** |
| DNBseq-G400RS 15x | 150 bp | 31 | 0.9015 | 2.8404 | 14.17× |
| HiSeqX 15x (controls) | 150 bp | 2,662 | 0.9159 | 2.9142 | 18.72× |
| NovaSeq 30x | 151 bp | 50 | 0.9179 | 2.9528 | 30.93× |

The 100 bp T7 library reaches a **higher** field depth than the 150 bp controls at
essentially the same measured depth. The informational worry — that 30 matching bases
discriminate less than 45 — does not show up at this depth. Consistent with that, the
fraction of reads surviving `-m 100` over the MHC window is 99.622 % for T7 against
99.628 % for HiSeqX: these CRAMs are fixed-length and untrimmed, so the threshold has no
left tail to bite into.

**Depth — the groups are aligned, so it cannot be the explanation.** Median measured depth
is 18.72× in controls and 18.65× in cases, Mann–Whitney *P* = 0.194. The `30x`/`15x` labels
are targets and disagree with the measurement; `Observed_Depth` is what is used everywhere
here.

**What is left, and it survives every control available.** Controls carry 4.63 % of called
chromosomes on an allele the reference panel never lists, against 3.53 % in cases — a ratio
of 1.31, Fisher *P* = 1.89e-07. Two things could explain that away and neither does:

- **Not depth.** Restricted to the 17–21× window where HiSeqX and DNBSeq-T7 both sit at a
  median near 18.5×, the shares are 4.48 % against 3.16 %, Fisher *P* = 3.58e-07.
- **Not the reference panel.** The same calls scored against the 105-individual 1000
  Genomes JPT panel give 7.42 % against 5.40 % — the *level* moves by more than half
  again, and the *ratio* moves from 1.31 to 1.37. The level is a property of the panel;
  the contrast is not.

Depth still matters where phenotype is held fixed. Within cases only, the unconfirmed share
falls with depth: DNBseq-G400RS 15x at 14.17× gives 10.66 %, DNBSeq-T7 at 18.43× gives
3.36 %, NovaSeq at 30.93× gives 1.32 %. The 31 DNBseq-G400RS 15x cases are the group to
watch — by a wide margin the worst-typed group in the study.
`05.qc/typing_confound.tsv` and `figures/typing_confound.png` carry all of it.

**The ancestry restriction did not create this and did not remove it.** Before the cohort
filter, over 3,569 samples, the same contrast was 4.74 % against 3.57 %, a ratio of 1.33.
Restricting to `full_mainland` removes the one explanation that would have been benign —
that the controls' unlisted alleles are simply not Japanese — and the contrast survives
essentially unchanged. That makes it a stronger finding, not a weaker one. The
3,569-sample tables are kept under `results/_superseded.3569/` and the ratios above are
recomputed from them by `verify.sh`, not carried over.

## Known limitations, in order of how much they should worry a reader

The deferred defects, with their evidence, are in [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md).
The one that changes a number is **§1**: DRB3/4/5 hemizygotes are encoded as homozygotes,
so those three genes' dosages are inflated. The list below is about the study rather than
the code.

1. **There is no PER-SAMPLE accuracy measurement.** The QC says whether a call was
   produced and how deeply it resolved, never whether a given sample's call is right.
   `05.qc/allele_frequency_check.tsv` closes half of this: the controls' 2-field
   frequencies are compared to a published Japanese panel, which would move if the typing
   were systematically mis-assigning alleles. It compares *population frequencies*, so it
   still cannot say a given sample is right. The measurement that would settle it — typing
   samples with a published type — is scoped and deferred in
   [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) §2.
2. **The typed set is the analysis cohort, so it cannot be re-used for a wider one.**
   Restricting typing to `full_mainland` gives every table one denominator, at the cost
   that a future analysis on a different sample set must re-derive its tables. The raw
   HLA-HD output for the 468 excluded samples is kept under `results/_superseded.3569/`
   precisely so that this costs a re-run of `collect_alleles.py` and not a re-run of
   HLA-HD.
3. **HLA-HD's licence is per-registrant and this uses the group install.** The binary
   cannot change silently, but for publication the component should run against a copy
   obtained under our own registration.
4. **Ambiguous calls are reported at the depth their candidates agree on**, which is
   sometimes 2 fields where a 3-field answer might exist. The full candidate list is kept,
   so nothing is lost, but a downstream analysis that assumes uniform resolution is wrong.
5. **14 of the 33 typed genes cannot carry residues** — IPD-IMGT/HLA publishes no protein
   alignment for them. Their allele calls are still in `allele_calls.tsv`.
6. **111 CRAMs had no index and were indexed by this run.** They live in a directory named
   `cram4temporary`; if it is cleared, those samples become unrunnable and the cohort
   silently shrinks. The manifest step fails loudly rather than proceeding.
