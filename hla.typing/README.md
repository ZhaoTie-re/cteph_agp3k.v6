# hla.typing

Read-based HLA allele typing and amino-acid residue calling for the **`full_mainland`
analysis cohort**. CRAM → HLA-region reads → HLA-HD allele calls → allele dosages and a
residue matrix at **IMGT numbering**.

→ Method and rationale: **[docs/METHODS.md](docs/METHODS.md)** ·
Output tree and data dictionary: **[docs/OUTPUTS.md](docs/OUTPUTS.md)** ·
This study's numbers: **[docs/STUDY_NOTES.md](docs/STUDY_NOTES.md)** ·
Known defects, deliberately left in: **[docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md)** ·
Delivery report: **[report/hla_typing_report.qmd](report/hla_typing_report.qmd)** ·
Shared visual grammar: **[../analysis/_shared/docs/FIGURES.md](../analysis/_shared/docs/FIGURES.md)**

## Which samples this covers

**Three filters, in order, and they are not the same kind of decision.**

| stage | rule | n | AGP3K | PH |
|---|---|---|---|---|
| the v6 cohort | `Flag_JHRPv6 = True` in the workbook | 3,592 | 3,135 | 457 |
| analysable | in `wgs.auto.par`'s `sample_qc.keep.id` | 3,569 | 3,117 | 452 |
| **comparable — the analysis cohort** | in `PopGMM_output/full_mainland.fid_iid.txt` | **3,101** | **2,662** | **439** |

The workbook drops rows that are not v6. Sample QC drops 23 heterozygosity outliers, which
for HLA typing are exactly the samples least likely to resolve into two clean haplotypes.
The third filter drops a further 468 and is an **ancestry** decision, not a quality one:
`full_mainland` is the union of 17 components of a Gaussian mixture fitted in BBJ
principal-component space, and `PopGMM_output/keep_list_summary.tsv` is its only in-repo
provenance.

**Why the cohort and not everything that typed.** The only external check available here
asks whether an allele we called appears in a mainland Japanese reference panel. A control
who is not mainland Japanese fails that check for a reason that has nothing to do with
typing. Typing the whole QC-passing set and subsetting downstream is the other valid
design; this one is chosen so that every table published here has one denominator.

The 468 samples the ancestry filter removes were typed by an earlier run and are **not
deleted**. Their reads, their HLA-HD output and the tables built from them are kept under
`results/_superseded.3569/`, and `verify.sh` §23 checks that they are all still there.

## What this component does *not* establish

**Read [docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md) before using the output.** Two
things there change how a result should be read:

- **§3** — **4.63 % of control chromosomes carry an allele the reference panel has never
  observed, against 3.53 % in cases**, a ratio of 1.31. It survives matching on measured
  depth, and it survives changing the reference panel: against a much smaller panel the
  level moves to 7.42 % and 5.40 % while the ratio stays at 1.37. So it is neither depth
  nor an artefact of what the panel happens to carry. What is left is platform-and-cohort,
  which is perfectly confounded with phenotype. A common allele depleted more in controls
  reads as enrichment in cases — a false risk association from the typing alone. Check rare
  rows of `allele_dosage.tsv` against `05.qc/allele_pgroup_map.tsv` (`in_reference = 0`)
  before testing them. Residue tests are affected least.
- **§1** — the DRB3, DRB4 and DRB5 columns of `allele_dosage.tsv` and `residue_dosage.tsv`
  encode hemizygotes as homozygotes, so those three genes' dosages are not trustworthy. The
  other 16 genes are unaffected.

Nothing about association. It produces genotypes; whether a residue is associated with
anything is for the component that reads them. It also does not establish typing
**accuracy** — there is no external truth set here, so the QC measures whether a call was
produced and how deeply it resolved, not whether it is right.

The platform confound this study sits on applies in full: cases were sequenced on
DNBSeq/NovaSeq, controls entirely on HiSeqX, with zero overlap. **A case/control
difference in typing quality cannot be attributed**, which is why every QC metric here is
reported per platform first.

## The contract with upstream

| what | from | value in this study |
|---|---|---|
| the v6 cohort | `info/…xlsx`, `Flag_JHRPv6 == True` | 3,592 |
| which of them is analysable | `wgs.auto.par/results/07_sample_qc/.../sample_qc.keep.id` | 3,569 |
| **which of those is comparable** | **`PopGMM_output/full_mainland.fid_iid.txt`** | **3,101** — 2,662 AGP3K, 439 PH |
| CRAM path | the same sheet, `Cram_Path` | 19–67 GB each, ~88 TB total |
| the CRAM index | **resolved by `realpath`**, not by appending `.crai` | most paths are symlinks; 111 CRAMs have no index and are indexed by this run |
| reference | `nagasaki_pipeline/data/hs38DH.fa` | 3,366 sequences, matching the CRAM `@SQ` exactly |
| allele dictionary | pinned `software/hlahd.dictionary/IMGT-HLA-3.64.0_N150` | IPD-IMGT/HLA 3.64.0, 46,005 alleles |
| residue alignments | pinned `software/hlahd.dictionary/IMGT-alignments-3.64.0` | the same release, 19 usable genes |
| gene list | pinned `software/hlahd.dictionary/HLA_gene.split.3.50.0.txt` | 33 genes; defines what every result is validated against |
| primary frequency reference | pinned `../jMorp_HLA_types/HLA_allele_frequencies_61K.txt` | ToMMo jMorp 61KJPN-HLA, 61,424 Japanese individuals, 13 loci |
| secondary frequency reference | pinned `../1KG_HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt` | 1000 Genomes JPT, 105 individuals, 5 loci |
| P-group definitions | pinned `../jMorp_HLA_types/hla_nom_p.txt` | IPD-IMGT/HLA 3.64.0, the allele-name translation the reference needs |

Every reference above is a **pinned copy**. The directories they came from belong to other
users or to a public FTP server and change without notice; a reference that moves under a
study is not a reference. `params.reference` is the one exception, and deliberately: it is
the hs38DH the CRAMs were *aligned to*, so it is fixed by the data rather than chosen here.

**Both frequency panels run on every execution.** That is not redundancy — see METHODS §10.

## Run it

```bash
source activate dsl2
cd hla.typing
nextflow run hla.typing.nf -resume
nextflow run hla.typing.nf --MaxSamples 50            # a pilot, spread over every platform
nextflow run hla.typing.nf --cohort_keep ''           # type every QC-passing sample instead
nextflow run hla.typing.nf --cohort_keep other.fid_iid.txt --CohortName narrow_mainland
```

Analysis scripts run under the `cteph_geno_pro` conda env, activated inside each process.
`bowtie2` and `samtools` come from `software/`; `hlahd.sh` from the group install.

**Run the pilot first.** The full cohort is ~480 CPU-hours, and the pilot is what confirms
the resource ladder and the QC gate before committing it.

**Run `./verify.sh` before and after.** Sections 1–10 check the source tree and need no
results; 11–23 check the tables. It is the only thing standing between this component and
the documentation rot it has had before — see METHODS §14.

### Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `sample_qc_keep` | `wgs.auto.par/…/sample_qc.keep.id` | the samples that passed sample QC; the component inherits this decision |
| `cohort_keep` | `PopGMM_output/full_mainland.fid_iid.txt` | the analysis cohort, applied after sample QC. `''` types every QC-passing sample |
| `CohortName` | `full_mainland` | recorded in `run_manifest.json`; does not affect the run |
| `MaxSamples` | `0` | 0 = all; a positive value takes a pilot spread across every platform |
| `MhcStart` / `MhcEnd` | 29,540,000 / 33,420,000 | the chr6 extraction window |
| `MinReadLength` | 100 | `hlahd.sh -m`; reads shorter than this are ignored |
| `CuttingRate` | 0.95 | `hlahd.sh -c` |
| `ResidueGenes` | 19 genes | those with an IMGT protein alignment |
| `MinResidueCount` | 1 | drop a residue or allele seen on fewer chromosomes |
| `AlleleFieldDepth` | 2 | resolution of `allele_dosage.tsv`, and of the frequency check |
| `RefLociJmorp` | 13 loci | the loci the primary panel carries |
| `RefLoci1kg` | 5 loci | the loci the secondary panel carries. **One list per panel** — sharing them is METHODS §10's measured failure |
| `ControlGroup` | `AGP3K` | the `group` value treated as a population sample |
| `maxForksExtract` | 48 | concurrent CRAM readers; the only concurrency limit |

## The fourteen steps

| # | process | writes to |
|---|---|---|
| 1 | `BUILD_MANIFEST` | `00.manifest/` — the sample set after all three filters, every CRAM and index resolved and verified before a job is submitted |
| 2 | `CONTIG_LIST` | `00.manifest/` — the 543 extraction regions, computed once for the run |
| 3 | `INDEX_CRAM` | `00.manifest/crai/` — an index for each CRAM that has none |
| 4 | `EXTRACT_READS` | `01.reads/` — CRAM to the FASTQ pair HLA-HD reads |
| 5 | `HLAHD` | `02.typing/` — typing **and the gate on it**, with the whole alignment tree kept |
| 6 | `COLLECT_ALLELES` | `03.alleles/` — allele / ambiguity / status tables |
| 7 | `BUILD_RESIDUE_REF` | `04.residues/` — IMGT alignments to a residue reference |
| 8 | `RESIDUE_MATRIX` | `04.residues/` — allele dosage, residue diplotype, residue dosage |
| 9 | `TYPING_QC` | `05.qc/` — per sample, platform, group and gene |
| 10 | `ALLELE_FREQ_CHECK` | `05.qc/` — **one task per reference panel**; control frequencies plus `allele_pgroup_map.tsv` for the association |
| 11 | `PLOT_TYPING_QC` | `figures/typing_qc.png` |
| 12 | `PLOT_ALLELE_FREQ` | `figures/allele_frequency.png` and `figures/allele_frequency_loci.png` |
| 13 | `PLOT_TYPING_CONFOUND` | `figures/typing_confound.png` and `05.qc/typing_confound.tsv` |
| 14 | `WRITE_RUN_MANIFEST` | `_run_info/` |

**Validation is the last line of step 5, and that placement is the point.** HLA-HD exits 0
whatever happens inside it, so a truncated result is a silent success. `validate_typing.py`
runs inside the `HLAHD` task, so a truncation is a non-zero exit **of the typing task** —
which is what the memory ladder retries, and which means nothing truncated is ever
published. As a separate downstream process its failure could only kill the run, after the
bad result had already been written. See METHODS §4.

### Concurrency

`EXTRACT_READS` and `INDEX_CRAM` are capped at `params.maxForksExtract`. They are the only
I/O-bound stages, and running all of them at once would pull ~88 TB off a **shared** Lustre
filesystem as fast as the scheduler can start them — a cost paid by everyone on it. `HLAHD`
is compute-bound and is left to the scheduler.

## Scripts

```
scripts/
  build_manifest.py      sample sheet + two keep lists -> verified CRAM + index
                         manifest; aborts at launch rather than dropping a sample
  contig_list.py         hs38DH.fa.fai -> the extraction regions
  extract_hla_reads.py   CRAM -> collated FASTQ pair, every exit code checked
  validate_typing.py     the gate: makes a silent HLA-HD failure non-zero. Runs
                         INSIDE the HLAHD task, never downstream of it
  collect_alleles.py     per-sample results -> allele / ambiguity / status tables
  build_residue_ref.py   IMGT *_prot.txt -> (allele x IMGT position) reference
  residue_matrix.py      calls x reference -> allele dosage, residue diplotype and
                         0/1/2 residue dosage, plus the allele -> reference match map
  typing_qc.py           call rate, resolution, ambiguity and failure per platform,
                         per case/control group and per gene
  allele_freq_check.py   control allele frequencies vs one reference panel, matched
                         on IPD-IMGT P groups. Run once per panel
  build_facts.py         every number this component's prose is allowed to quote,
                         derived from the published tables into _run_info/facts.json
  plot_typing_qc.py      the COMPLETENESS figure
  plot_allele_freq.py    the frequency check: a summary figure and a per-locus record
  plot_typing_confound.py  the CONFOUNDING figure, under both reference panels
```

## Checks

`./verify.sh` is the check. It replaces the ad-hoc greps that used to live here, and it is
the only reason the numbers in these documents can be trusted: section 9 fails if any
number in any `.md` is neither in `_run_info/facts.json` nor in
`docs/NUMBERS_ALLOWED.md`, and section 8 fails if any figure caption states a number
instead of interpolating one.

```bash
./verify.sh                       # everything; 1-10 need no results
column -t results/05.qc/allele_frequency_summary.tsv
du -sh  results/                  # copied only
du -shL results/                  # following the symlinks into work/
```
