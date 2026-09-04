# STUDY NOTES — hla.typing

This study's configuration and its measured numbers. Method and rationale are in
[METHODS.md](METHODS.md).

## Configuration

| | |
|---|---|
| samples | **3,569** — `Flag_JHRPv6 == True` (3,592) minus the **23** heterozygosity outliers `wgs.auto.par` sample QC excluded |
| platforms | HiSeqX 15x 3,117 · DNBSeq-T7 30x 327 · NovaSeq 30x 51 · DNBSeq-G400RS 30x 43 · DNBseq-G400RS 15x 31 |
| excluded by sample QC | 23 — 18 AGP3K (HiSeqX), 5 PH (T7 4, NovaSeq 1); all `n_fail_het` |
| reference | `nagasaki_pipeline/data/hs38DH.fa` — 3,366 sequences, 525 `HLA*` contigs |
| extraction regions | **543** = MHC window + 16 `chr6_*_alt` + 525 `HLA*` + unmapped |
| typing | HLA-HD 1.7.1, `-t 4 -m 100 -c 0.95` |
| dictionary | pinned IPD-IMGT/HLA **3.64.0**, 46,005 alleles, `_N150` |
| residues | pinned IMGT protein alignments, **same release**, 19 genes |
| gene list | pinned `HLA_gene.split.3.50.0.txt`, 33 genes |
| frequency reference | pinned **jMorp 61KJPN-HLA** — 61,424 Japanese individuals, 122,848 chromosomes per locus, 13 loci, matched on IPD-IMGT P groups. The 1000 Genomes JPT panel (105 samples, 5 loci) is kept as an independent second reference |

## Sanity anchors — a refactor must not move these

The residue numbering is IMGT's, and these are how that is known:

| | |
|---|---|
| HLA-A positions 1–8 | **`GSHSMRYF`** — the canonical mature N-terminus |
| *DRB1* `*04:05:01` at 11/13/71/74 | **V / H / R / A** |
| *DRB1* `*15:01:01` at 11/13/71/74 | **P / R / A / A** |
| *DRB1* `*03:01:01` at 11/13/71/74 | **S / S / K / R** |
| `NXNE` in `residue_reference/` | **0 rows** — the signature of a fabricated residue |
| JPT reference frequencies | `A*24:02` 0.390 · `B*52:01` 0.124 · `C*07:02` 0.152 · `DRB1*09:01` 0.143 · `DQB1*06:01` 0.188 |
| the frequency check against itself | feeding the JPT panel in as if it were the calls gives `r = ρ = 1.000` and `max\|diff\| = 0.000` at all five loci |

Reference sizes: A 9,175 alleles × 390 positions (−30..360) · B 11,103 × 483 ·
DRB1 3,977 × 290 (−30..236) · DQB1 3,068 × 274 (−40..234).

## Measured on this data

**Extraction.** One HiSeqX control: 628,405 read pairs over the 543 regions, R1 and R2
counts equal, 6,731 singletons captured to their own file.

**Typing.** The same sample end to end: **8 min 14 s** wall clock at `-t 4`, peak RSS
**1.51 GB**, 33 genes typed. The reference implementation at `-t 1` takes 30–50 minutes.
Extrapolated: ~480 CPU-hours for the cohort, about 10 hours at concurrency 50.

**Output size, and why none of it is copied.** Per sample: ~90 MB of extracted FASTQ and
~1.2 GB of HLA-HD output, of which 1,174 MB is intermediate SAM (`mapfile/` 867 M, `intron/`
188 M, `exon/` 136 M). Across the cohort that is **321 GB + 4.2 TB**. All of it is kept and
all of it is **symlinked** into `work/`; only the tables, figures and run record are copied,
and those are under 100 MB. Copying the rest would duplicate 4.5 TB; deleting it would mean
re-reading 88 TB of CRAM to get back. The cost of the choice is that **`work/` must not be
cleaned** — see [OUTPUTS.md](OUTPUTS.md).

**Input.** CRAMs are 19–67 GB each, ~88 TB in total to read. Extraction is I/O bound.

## The confounding this component sits on

Platform is **perfectly separated from phenotype**: every control is HiSeqX, every case
is DNBSeq or NovaSeq. Nothing in this component can fix that, and no covariate can.

The two axes that could plausibly bias typing were checked before the component was built:

**Read length — measured, and the answer is no effect.** DNBSeq-T7 is a **100 bp** library
(327 of 452 cases) against 150 bp everywhere else. On the **full 3,569-sample run**:

| platform | read length | *n* | call rate | field depth | measured depth |
|---|---|---|---|---|---|
| DNBSeq-G400RS 30x | 150 bp | 43 | 0.9204 | 2.9255 | 36.4× |
| **DNBSeq-T7 30x** | **100 bp** | **327** | **0.9182** | **2.9310** | **18.4×** |
| DNBseq-G400RS 15x | 150 bp | 31 | 0.9015 | 2.8404 | 14.2× |
| HiSeqX 15x (controls) | 150 bp | 3,117 | 0.9162 | 2.9135 | 19.0× |
| NovaSeq 30x | 151 bp | 51 | 0.9174 | 2.9514 | 30.9× |

The 100 bp T7 library reaches a **higher** call rate and field depth than the 150 bp controls
at essentially the same measured depth. The informational worry (30 matching bases discriminate
less than 45) does not show up at this depth.

**Call rate and field depth are nearly saturated and nearly flat**, which is why they are weak
quality metrics here: call rate takes **five distinct values** across 3,569 samples (it is
*k*/19) and 97 % of samples sit on just two of them. What actually varies between platforms is
accuracy, and only the frequency check of §3 sees it.

Consistent with that, the fraction of reads surviving `-m 100` over `chr6:29.54–33.42 Mb` is
T7 **99.622 %** against HiSeqX **99.628 %** — these CRAMs are fixed-length and untrimmed, so
the threshold has no left tail to bite into.

**Depth — this is what actually correlates.** Field depth against `Observed_Depth`,
Spearman **ρ = +0.299, *P* = 0.035**; call rate ρ = +0.154, *P* = 0.286. The weakest platform
in the pilot is DNBseq-G400RS **15x** — the shallowest at 12.0× MHC coverage, and 31 cases.

The `30x`/`15x` labels are targets. Measured `Observed_Depth`: cases median **18.69×**,
controls **19.03×** — the main comparison is depth-aligned, because the T7 correction in the
workbook's edit history (a 100 bp library whose depth was computed as if reads were 150 bp)
brought 331 cases down to ~18.5×. Directly measured MHC coverage: G400RS-30x 41.5×/32.1× ·
NovaSeq 30.9× · HiSeqX 24.5×/22.7× · T7 19.9×/21.1× · G400RS-15x 12.0×. T7 carries about
**20 % more reads** than the controls despite lower coverage, because HLA-HD counts reads and
T7's are shorter.

**This depth-aligned reading was right, and OPEN_QUESTIONS §3 originally contradicted it.**
That section attributed the accuracy gap to depth by reading the `15x`/`30x` labels as
measurements; it has been corrected. Restricted to the 17–21× window where HiSeqX and T7 both
sit at a median of 18.6×, the unconfirmed share is **4.71 % against 3.13 %**, Fisher
*P* = 9.6 × 10⁻¹⁰ — the gap is **platform and cohort, at matched depth**, and it is not
separable from phenotype in these data.

Depth still matters where phenotype is held fixed. Within cases only: G400RS 14.2× → **10.7 %**
unconfirmed, T7 18.4× → 3.4 %, NovaSeq 30.9× → 1.5 %. The 31 DNBseq-G400RS 15x cases are the
group to watch — 7 % of cases, and by a wide margin the worst-typed group in the study.
`05.qc/typing_confound.tsv` and `figures/typing_confound.png` carry the numbers.

All of these fell on 2026-08-26 when the reference changed from 1000 Genomes JPT (105 people,
5 loci) to jMorp 61KJPN (61,424 people, 13 loci). **The typing did not change and neither did
the conclusion** — a better reference simply recognises more of what we called, so less of the
tail is a reference gap. OPEN_QUESTIONS §3 records the before-and-after.

## Known limitations, in order of how much they should worry a reader

The deferred defects, with their evidence, are in [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md). The
one that changes a number is **§1**: DRB3/4/5 hemizygotes are encoded as homozygotes, so those
three genes' dosages are inflated. The list below is about the study rather than the code.

1. **There is no PER-SAMPLE accuracy measurement.** The QC says whether a call was
   produced and how deeply it resolved, never whether a given sample's call is right.
   `05.qc/allele_frequency_check.tsv` closes half of this: the controls' 2-field
   frequencies at A/B/C/DQB1/DRB1 are compared to the published 1000 Genomes **JPT**
   panel, which would move if the typing were systematically mis-assigning alleles. It
   compares *population frequencies*, so it still cannot say a given sample is right, and
   it cannot reach DPB1, DQA1 or DRB3/4/5 at all — the loci typed least reliably here.
   The measurement that would settle it (typing the 104 JPT samples that have both a
   published type and a downloadable CRAM, ~1.57 TB) is scoped and deferred in
   [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) §2, together with the read-length experiment.
2. **HLA-HD's licence is per-registrant and this uses the group install.** The binary
   cannot change silently, but for publication the component should run against a copy
   obtained under our own registration.
3. **Ambiguous calls are reported at the depth their candidates agree on**, which is
   sometimes 2 fields where a 3-field answer might exist. The full candidate list is kept,
   so nothing is lost, but a downstream analysis that assumes uniform resolution is wrong.
4. **14 of the 33 typed genes cannot carry residues** — IPD-IMGT/HLA publishes no protein
   alignment for them. Their allele calls are still in `allele_calls.tsv`.
5. **111 CRAMs had no index and were indexed by this run.** They live in a directory named
   `cram4temporary`; if it is cleared, those samples become unrunnable and the cohort
   silently shrinks by 3 %. The manifest step fails loudly rather than proceeding.
