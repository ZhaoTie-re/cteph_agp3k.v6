# hla.typing

Read-based HLA allele typing and amino-acid residue calling, for every sample that is
**both** marked `Flag_JHRPv6 = True` in the workbook **and** passed `wgs.auto.par`'s
sample QC. CRAM → HLA-region reads → HLA-HD allele calls → allele dosages and a
residue matrix at **IMGT numbering**.

→ Method and rationale: **[docs/METHODS.md](docs/METHODS.md)** ·
Output tree and data dictionary: **[docs/OUTPUTS.md](docs/OUTPUTS.md)** ·
This study's numbers: **[docs/STUDY_NOTES.md](docs/STUDY_NOTES.md)** ·
Known defects, deliberately left in: **[docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md)** ·
Shared visual grammar: **[../analysis/_shared/docs/FIGURES.md](../analysis/_shared/docs/FIGURES.md)**

## What this component does *not* establish

**Read [docs/OPEN_QUESTIONS.md](docs/OPEN_QUESTIONS.md) before using the output.** Two
things there change how a result should be read:

- **§3** — **4.6 % of typed chromosomes carry an allele that a 61,424-person Japanese panel
  has never observed**, and that mass comes out of the common alleles (`A*24:02P` at 0.315
  here against a published 0.339). It **is** differential: controls sit at 4.74 % against
  cases at 3.57 %, and the gap survives matching on measured depth (4.71 % vs 3.13 %, Fisher
  *P* = 9.6 × 10⁻¹⁰), so it tracks platform-and-cohort, which is perfectly confounded with
  phenotype. A common allele depleted more in controls reads as enrichment in cases — a false
  risk association from the typing alone. Check rare rows of `allele_dosage.tsv` against
  `05.qc/allele_pgroup_map.tsv` (`in_reference = 0`) before testing them. Residue tests are
  affected least.
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
| **which of them is analysable** | **`wgs.auto.par/results/07_sample_qc/.../sample_qc.keep.id`** | **3,569** — 3,117 AGP3K, 452 PH; 23 heterozygosity outliers excluded |
| CRAM path | the same sheet, `Cram_Path` | 19–67 GB each, ~88 TB total |
| the CRAM index | **resolved by `realpath`**, not by appending `.crai` | 3,117 paths are symlinks; 111 CRAMs have no index and are indexed by this run |
| reference | `nagasaki_pipeline/data/hs38DH.fa` | 3,366 sequences, matching the CRAM `@SQ` exactly |
| allele dictionary | pinned `software/hlahd.dictionary/IMGT-HLA-3.64.0_N150` | IPD-IMGT/HLA 3.64.0, 46,005 alleles |
| residue alignments | pinned `software/hlahd.dictionary/IMGT-alignments-3.64.0` | the same release, 19 usable genes |
| gene list | pinned `software/hlahd.dictionary/HLA_gene.split.3.50.0.txt` | 33 genes; defines what every result is validated against |
| frequency reference | pinned `../jMorp_HLA_types/HLA_allele_frequencies_61K.txt` | ToMMo jMorp 61KJPN-HLA, 61,424 Japanese individuals, 13 loci |
| P-group definitions | pinned `../jMorp_HLA_types/hla_nom_p.txt` | IPD-IMGT/HLA 3.64.0, the allele-name translation the reference needs |

All four are **pinned copies**. The directories they came from belong to other users or
to a public FTP server and change without notice; a reference that moves under a study
is not a reference. `params.reference` is the one exception, and deliberately: it is the
hs38DH the CRAMs were *aligned to*, so it is fixed by the data rather than chosen here.

## Run it

```bash
source activate dsl2
cd hla.typing
nextflow run hla.typing.nf -resume
nextflow run hla.typing.nf --MaxSamples 50        # a pilot, spread over every platform
nextflow run hla.typing.nf --sample_qc_keep other.keep.id   # a different QC decision
nextflow run hla.typing.nf --MinReadLength 90     # override an HLA-HD knob
```

Analysis scripts run under the `cteph_geno_pro` conda env, activated inside each process.
`bowtie2` and `samtools` come from `software/`; `hlahd.sh` from the group install.

**Run the pilot first.** The full cohort is ~480 CPU-hours, and the pilot is what confirms
the resource ladder and the QC gate before committing it.

### Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `sample_qc_keep` | `wgs.auto.par/…/sample_qc.keep.id` | the samples that passed sample QC; the component inherits this decision rather than making its own |
| `MaxSamples` | `0` | 0 = all; a positive value takes a pilot spread across every platform |
| `MhcStart` / `MhcEnd` | 29,540,000 / 33,420,000 | the chr6 extraction window |
| `MinReadLength` | 100 | `hlahd.sh -m`; reads shorter than this are ignored |
| `CuttingRate` | 0.95 | `hlahd.sh -c` |
| `ResidueGenes` | 19 genes | those with an IMGT protein alignment |
| `MinResidueCount` | 1 | drop a residue or allele seen on fewer chromosomes |
| `AlleleFieldDepth` | 2 | resolution of `allele_dosage.tsv`, and of the frequency check |
| `TruthPanel` | jMorp 61KJPN-HLA | the reference the frequencies are compared to |
| `TruthFormat` | `jmorp_long` | `1kg_wide` switches back to the 1000 Genomes panel |
| `TruthPopulation` | `JPT` | `1kg_wide` only: the population within that panel |
| `ControlGroup` | `AGP3K` | the `group` value treated as a population sample |
| `maxForksExtract` | 48 | concurrent CRAM readers; the only concurrency limit |

## The fourteen steps

| # | process | writes to |
|---|---|---|
| 1 | `BUILD_MANIFEST` | `00.manifest/` — the sample set, every CRAM and index resolved and verified before a job is submitted |
| 2 | `CONTIG_LIST` | `00.manifest/` — the 543 extraction regions, computed once for the run |
| 3 | `INDEX_CRAM` | `00.manifest/crai/` — an index for each of the 111 CRAMs that has none |
| 4 | `EXTRACT_READS` | `01.reads/` — CRAM to the FASTQ pair HLA-HD reads |
| 5 | `HLAHD` | `02.typing/` — typing **and the gate on it**, with the whole alignment tree kept |
| 6 | `COLLECT_ALLELES` | `03.alleles/` — allele / ambiguity / status tables |
| 7 | `BUILD_RESIDUE_REF` | `04.residues/` — IMGT alignments to a residue reference |
| 8 | `RESIDUE_MATRIX` | `04.residues/` — allele dosage, residue diplotype, residue dosage |
| 9 | `TYPING_QC` | `05.qc/` — per sample, platform, group and gene |
| 10 | `ALLELE_FREQ_CHECK` | `05.qc/` — control frequencies against jMorp 61KJPN-HLA, plus `allele_pgroup_map.tsv` for the association |
| 11 | `PLOT_TYPING_QC` | `figures/typing_qc.png` |
| 12 | `PLOT_ALLELE_FREQ` | `figures/allele_frequency.png` |
| 13 | `PLOT_TYPING_CONFOUND` | `figures/typing_confound.png` |
| 14 | `WRITE_RUN_MANIFEST` | `_run_info/` |

**Validation is the last line of step 5, and that placement is the point.** HLA-HD exits 0
whatever happens inside it, so a truncated result is a silent success. `validate_typing.py`
runs inside the `HLAHD` task, so a truncation is a non-zero exit **of the typing task** —
which is what the memory ladder retries, and which means nothing truncated is ever
published. As a separate downstream process its failure could only kill the run, after the
bad result had already been written. See METHODS §4.

### Concurrency

`EXTRACT_READS` and `INDEX_CRAM` are capped at `params.maxForksExtract` (48). They are the
only I/O-bound stages, and 3,569 of them at once would pull ~88 TB off a **shared** Lustre
filesystem as fast as the scheduler can start them — a cost paid by everyone on it. `HLAHD`
is compute-bound and is left to the scheduler.

## Scripts

```
scripts/
  build_manifest.py      sample sheet -> verified CRAM + index manifest; aborts at launch
  contig_list.py         hs38DH.fa.fai -> the 543 extraction regions
  extract_hla_reads.py   CRAM -> collated FASTQ pair, every exit code checked
  validate_typing.py     the gate: makes a silent HLA-HD failure non-zero. Runs
                         INSIDE the HLAHD task, never downstream of it
  collect_alleles.py     per-sample results -> allele / ambiguity / status tables
  build_residue_ref.py   IMGT *_prot.txt -> (allele x IMGT position) reference
  residue_matrix.py      calls x reference -> allele dosage, residue diplotype and
                         0/1/2 residue dosage, plus the allele -> reference match map
  typing_qc.py           call rate, resolution, ambiguity and failure per platform,
                         per case/control group and per gene
  allele_freq_check.py   control allele frequencies vs jMorp 61KJPN-HLA (61,424
                         Japanese), matched on IPD-IMGT P groups. Also writes the
                         allele -> P group -> reference frequency crosswalk
  plot_typing_qc.py      the COMPLETENESS figure: field resolution, ambiguity and
                         the four per-gene outcomes, each as a composition
  plot_typing_confound.py  the CONFOUNDING figure: unconfirmed-allele share against
                         measured depth and platform, incl. the depth-matched window
  plot_allele_freq.py    the frequency-check figure
```

## Checks worth running after a run

```bash
# nothing failed
awk -F'\t' 'NR>1 && $4!="COMPLETED" && $4!="CACHED"' results/_run_info/trace.txt

# no silently truncated typing — nothing truncated should ever have been published
grep -rl "Couldn't read result file" results/02.typing/ || echo "clean"

# every sample has all 33 genes
awk 'END{print NR" gene rows in the last result file"}' \
  results/02.typing/*/result/*_final.result.txt | tail -1

# the residue reference is not fabricated: the NXNE signature of a failed
# blastx round-trip must not appear
grep -rc $'N\tX\tN\tE' results/04.residues/residue_reference/ | grep -v ':0' || echo "clean"

# no allele failed to match the reference at any depth (the map is the successes)
awk 'END{print NR-1" allele(s) with no match at all"}' results/04.residues/allele_unmatched.tsv

# the external check: control frequencies against jMorp 61KJPN-HLA
column -t results/05.qc/allele_frequency_summary.tsv

# what is copied stays small; what is symlinked is large and lives in work/
du -sh  results/       # copied only
du -shL results/       # following the symlinks
```
