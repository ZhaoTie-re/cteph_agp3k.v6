# Outputs — hla.typing

Every path below is under `results/`.

```
00.manifest/                   copied
  sample_manifest.tsv          the sample set, with every CRAM and index resolved
  contig_list.txt              the 543 extraction regions
  crai/<sample>.cram.crai      indexes built by this run for CRAMs that had none
01.reads/                      SYMLINKED   ~90 MB/sample, 321 GB total
  <sample>.hla.R{1,2}.fastq.gz         the reads HLA-HD typed
  <sample>.hla.singleton.fastq.gz      reads whose mate did not survive extraction
02.typing/<sample>/            SYMLINKED   ~1.2 GB/sample, 4.2 TB total
  result/                      the calls — <s>_final.result.txt, <s>_<GENE>.est.txt
  log/                         <s>_<GENE>.log, one per gene
  hlahd.run.log                the run log; the validator greps it for kill lines
  typing_validation.tsv        what the validator saw for this sample
  mapfile/ exon/ intron/ maplist/      the bowtie2 alignments the calls came from
03.alleles/                    copied
  allele_calls.tsv             one row per sample, three columns per gene
  allele_ambiguity.tsv         every ambiguous call, with all its candidates
  typing_status.tsv            per-sample counts of each outcome
04.residues/
  residue_reference/<gene>.tsv (allele x IMGT position) — SYMLINKED, 34 MB of
                               derived reference, the same for a given IMGT release
  allele_dosage.tsv            0/1/2 per 2-field allele — the classical allele test
  residue_diplotype.tsv · residue_dosage.tsv · residue_sites.tsv
  residue_reference_summary.tsv
  allele_match_map.tsv         how each allele reached the reference, once per allele
  allele_unmatched.tsv         only the alleles that matched NOTHING       all copied
05.qc/                         copied
  typing_qc_sample.tsv · typing_qc_platform.tsv · typing_qc_group.tsv
  typing_qc_gene.tsv
  allele_frequency_check.tsv · allele_frequency_summary.tsv
figures/                       copied
  typing_qc.png · typing_qc.md
  allele_frequency.png · allele_frequency.md
_run_info/                     copied
  trace.txt · report.html · timeline.html · dag.html · run_manifest.json
```

## Copied, symlinked, and the one thing that follows from it

**Everything a reader opens is copied**: the tables, the figures, the manifest, the run
record. Together they are under 100 MB and they survive anything.

**Everything large is symlinked and nothing is thrown away**: 321 GB of extracted reads,
4.2 TB of alignments, 34 MB of residue reference. Copying them would cost 4.5 TB to duplicate
data that already exists in `work/`; deleting them would mean re-reading 88 TB of CRAM to get
back to the same place. The symlink is the only option that costs nothing and loses nothing.

**Extraction is fork-limited.** `EXTRACT_READS` and `INDEX_CRAM` run at most
`params.maxForksExtract` (48) at a time. They are the only I/O-bound stages, and 3,569 of
them at once would pull ~88 TB off a shared filesystem as fast as the scheduler allows — a
cost paid by everyone on it, not only by this run.

**So: do not clean `work/`.** Cleaning it leaves every table and figure intact and empties
`01.reads/`, `02.typing/` and `04.residues/residue_reference/`. That is survivable — the
study's conclusions are all in the copied files — but re-typing with a newer IMGT release, or
validating against a second tool, would start again from the CRAMs.

## `sample_manifest.tsv`

| column | meaning |
|---|---|
| `sample_id` | `ID_JHRPv6` |
| `cram` | `Cram_Path` from the workbook, as written |
| `crai` | the index, resolved through `realpath`. **Empty means no index exists anywhere**, and the workflow builds one into `00.manifest/crai/` |
| `platform`, `group`, `sex`, `observed_depth` | carried from the workbook for QC stratification |

## `allele_calls.tsv`

One row per sample. Three columns per gene:

| column | meaning |
|---|---|
| `<GENE>_1`, `<GENE>_2` | the two allele calls, `HLA-` stripped. Empty when the state is `not_typed` or `failed` |
| `<GENE>_state` | `called` · `hemizygous` · `not_typed` · `failed` — see METHODS §8. `failed` is a pipeline failure, not biology |

## `allele_ambiguity.tsv`

| column | meaning |
|---|---|
| `n_candidates` | how many alleles HLA-HD listed (always an even number of pair members) |
| `candidates` | all of them, `;`-joined, in HLA-HD's order |
| `reported_1`, `reported_2` | what `allele_calls.tsv` carries |
| `agreed_field_depth` | the depth at which every candidate pair is the same genotype; 0 means they disagree at the first field |

## `residue_reference/<gene>.tsv`

Rows are alleles, columns are **IMGT positions**. Position 1 is the first residue of the
mature protein; the leader peptide runs negative and there is no zero. A cell is a
one-letter residue, `.` where the allele has no sequence there, or `*` where IMGT records
it as unknown.

## `allele_dosage.tsv`

One row per sample, one column per **2-field allele**, value 0/1/2 — the classical HLA
allele test, and the level most papers report first. Column names are the alleles
themselves (`A*24:02`, `DRB1*15:02`), so the gene is in the name.

Covers **every gene that produced calls**, not only the 19 with a residue reference. An
allele seen on fewer than `MinResidueCount` chromosomes is dropped, as for residues. Empty
where the call was `failed` or `not_typed`.

A gene's columns sum to `2 × (called + hemizygous)` from `typing_status.tsv`, because a
hemizygous call is written as the same allele twice. For A/B/C/DQA1/DQB1/DP that is what a
homozygote is; **for DRB3, DRB4 and DRB5 it is wrong** in the same way it is wrong in
`residue_dosage.tsv` — see [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) §1.

## `allele_match_map.tsv` and `allele_unmatched.tsv`

| column | meaning |
|---|---|
| `gene`, `allele` | the call, as HLA-HD wrote it |
| `matched_to` | the reference allele its residues were read from |
| `how` | `exact`, `prefix(N)` — N reference alleles share this prefix, the alphabetically first was taken — or `truncated_to_D` |
| `n_chromosomes` | how many chromosomes across the cohort carried this call |

**One row per allele, not per sample.** How a call matches the reference is a property of
the allele; recording it per *(sample, allele)* produced ~120,000 rows for this cohort of
which essentially all were ordinary successful prefix matches.

`allele_unmatched.tsv` holds only the calls that matched **nothing at any depth**, one row
per *(sample, gene, allele)* so the sample can be traced. It is normally empty apart from
its header, and a row in it means a genotype was lost.

## `residue_diplotype.tsv` and `residue_dosage.tsv`

| | rows | columns | value |
|---|---|---|---|
| `residue_diplotype.tsv` | samples | `<GENE>:<position>` | the two residues, **sorted**, e.g. `FH` |
| `residue_dosage.tsv` | samples | `<GENE>:<position>:<residue>` | 0, 1 or 2 |

**The DRB3, DRB4 and DRB5 columns of the dosage matrix are known to be wrong.** A hemizygote —
one copy of the gene, which is the common case for these three haplotype-dependent paralogues
— is encoded as a homozygote, so its dosages read 2 where they should read 1. See
[OPEN_QUESTIONS.md §1](OPEN_QUESTIONS.md). The allele calls are unaffected, and so are the
other 16 genes.

Both carry only positions polymorphic in this cohort. In the dosage matrix the columns for
one position sum to 2 in every sample with a call, which is worth asserting after a run.

## `residue_sites.tsv`

| column | meaning |
|---|---|
| `positions` | IMGT positions in the reference for this gene |
| `polymorphic` | how many of them vary in this cohort — the testable sites the gene contributes |
| `note` | `no IMGT protein alignment` for the 14 genes that cannot carry residues |

## `allele_frequency_check.tsv` and `allele_frequency_summary.tsv`

The external check: control allele frequencies against the published 1000 Genomes **JPT**
panel, at A, B, C, DQB1 and DRB1 — the five loci the panel carries, at two fields.

`allele_frequency_check.tsv`, one row per *(gene, allele)* seen in either source:

| column | meaning |
|---|---|
| `count_ctrl`, `n_chr_ctrl`, `freq_ctrl` | observed in the control group (`params.ControlGroup`) |
| `ci_lo_ctrl`, `ci_hi_ctrl` | 95 % Wilson interval — a normal approximation goes below zero for the rare alleles that are most of an HLA table |
| `count_case`, `n_chr_case`, `freq_case` | the cases, reported beside and **not** compared |
| `count_ref`, `n_chr_ref`, `freq_ref` | the JPT reference |
| `diff_ctrl_ref` | `freq_ctrl − freq_ref` |
| `fisher_p_ctrl_ref` | Fisher exact against the reference **counts**, unadjusted. The reference is itself an estimate from 105 samples; treating it as known would make every difference look significant |
| `in_ref_only`, `in_obs_only` | present in one source and absent from the other |

`allele_frequency_summary.tsv`, one row per locus: `n_alleles`, both denominators,
`pearson_r` and `spearman_rho` between the frequency vectors, `max_abs_diff` and the allele
where it occurs, and the two directional counts. **`r` is dominated by the common alleles
and `ρ` weights every allele equally**, so the two disagreeing localises the disagreement
to the tail.

This compares **population frequencies**, not genotypes. It cannot say a given sample was
typed correctly. See [METHODS.md](METHODS.md) §11 and
[OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) §2.

## `typing_qc_platform.tsv` and `typing_qc_group.tsv`

The same metric set, aggregated two ways. **`typing_qc_platform.tsv` is the table to
read.** Platform is fully confounded with phenotype in this study — cases on DNBSeq and
NovaSeq, controls entirely on HiSeqX, zero overlap — so platform is the only axis on which
a technical gradient is visible as itself.

`typing_qc_group.tsv` gives the case/control split because a reader will ask for it. It is
a **description, not a test**: with no platform shared between the groups, a difference in
it has no attributable cause.

## `run_manifest.json`

What actually ran: the sample sheet, the sample-QC keep list and `MaxSamples`; the
reference and MHC window; the HLA-HD binary and the **pinned** dictionary, gene-split and
IMGT alignment paths; `-m` and `-c`; the residue gene list and minimum count; the allele
field depth; and the pinned frequency reference with its population and control group.
