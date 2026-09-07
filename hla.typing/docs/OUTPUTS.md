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
  allele_frequency_check.tsv    one row per (locus, allele): ours, cases, reference
  allele_frequency_summary.tsv  one row per locus: r, rho, max|diff|
  allele_frequency_sample.tsv   one row per sample: chromosomes called, unconfirmed
  allele_pgroup_map.tsv         allele -> P group -> in_reference; the association
                                reads this one
  *.1kg.tsv                     the same four against the SECONDARY panel. Same
                                schema, different reference and locus list
  typing_confound.tsv           written by PLOT_TYPING_CONFOUND: every stratum
                                behind figures/typing_confound.png, both panels
figures/                       copied
  typing_qc.png · typing_qc.md
  allele_frequency.png · allele_frequency.md            the summary
  allele_frequency_loci.png · allele_frequency_loci.md  the per-locus record
  typing_confound.png · typing_confound.md
_run_info/                     copied
  trace.txt · report.html · timeline.html · dag.html · run_manifest.json
  facts.json                   every number this component's prose may quote,
                               derived from the tables above by build_facts.py --
                               plus the jMorp allele list itself, which is the only
                               way to measure what NOT harmonising would have cost.
                               verify.sh section 9 fails on any number in a .md
                               that is not in here or in docs/NUMBERS_ALLOWED.md
_superseded.3569/              NOT AN OUTPUT — the record of what the cohort filter
                               removed. The 468 samples' reads and HLA-HD output
                               (moved, still symlinks into work/) and the tables
                               built before the filter existed, with their md5s.
                               verify.sh section 23 checks it is intact.
```

## Copied, symlinked, and the one thing that follows from it

**Everything a reader opens is copied**: the tables, the figures, the manifest, the run
record. Together they are under 100 MB and they survive anything.

**Everything large is symlinked and nothing is thrown away**: 321 GB of extracted reads,
4.2 TB of alignments, 34 MB of residue reference. Copying them would cost 4.5 TB to duplicate
data that already exists in `work/`; deleting them would mean re-reading 88 TB of CRAM to get
back to the same place. The symlink is the only option that costs nothing and loses nothing.

**Extraction is fork-limited.** `EXTRACT_READS` and `INDEX_CRAM` run at most
`params.maxForksExtract` (48) at a time. They are the only I/O-bound stages, and all of
them at once would pull ~88 TB off a shared filesystem as fast as the scheduler allows — a
cost paid by everyone on it, not only by this run.

**So: do not clean `work/`.** Cleaning it leaves every table and figure intact and empties
`01.reads/`, `02.typing/` and `04.residues/residue_reference/`. That is survivable — the
study's conclusions are all in the copied files — but re-typing with a newer IMGT release, or
validating against a second tool, would start again from the CRAMs.

## READ `sample_id` AS A STRING, ALWAYS

**66 of the 3,101 sample ids carry leading zeros** — `0000000123`, not `123`. Every file in
this tree preserves them, but `pandas.read_csv` without `dtype` coerces the column to an integer
and silently drops the zeros. A join keyed on it then loses those 66 samples with no error and
no warning: an inner merge simply returns 3,035 rows.

```python
pd.read_csv(path, sep='\t', dtype={'sample_id': str})          # tables
pd.read_csv(path, sep='\t', dtype={'sample_id': str}).set_index('sample_id')   # matrices
```

This bit the audit of 2026-08-26 before it bit anyone else: two independently computed allele
counts appeared to disagree on 201 of 1,543 alleles, and the disagreement was entirely this.

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

Both carry only positions polymorphic in this cohort.

**The columns for one position do NOT always sum to 2.** They sum to 2 for 98.61 % of
(sample, position) pairs, to 1 for 1.4 % and to 0 for 0.16 %. IMGT's protein alignment writes
`*` for an unsequenced residue and `.` for an alignment gap, and an allele carrying either
contributes no residue at that position — HLA-A position -22 is `*` for 2,703 of the alignment's
alleles. So a residue dosage is the count among chromosomes with a DETERMINED residue there, not
among all typed chromosomes.

Over the 1,055 positions at the 16 genes with usable dosages the median undetermined rate is
0.0322 %, but **143 positions exceed 5 %** and DQA1:56 reaches 23.73 %. `assoc_hla` carries this as
a per-position QC column and holds low-determination positions out of its omnibus test.

## `residue_sites.tsv`

| column | meaning |
|---|---|
| `positions` | IMGT positions in the reference for this gene |
| `polymorphic` | how many of them vary in this cohort — the testable sites the gene contributes |
| `note` | `no IMGT protein alignment` for the 14 genes that cannot carry residues |

## `allele_frequency_check.tsv` and `allele_frequency_summary.tsv`

The external check: control allele frequencies against **jMorp 61KJPN-HLA** — 61,424 Japanese
individuals over 13 loci. Alleles are matched on IPD-IMGT P groups, which is how the reference
names them; see [METHODS.md](METHODS.md) §11.

`allele_frequency_check.tsv`, one row per *(gene, allele)* seen in either source:

| column | meaning |
|---|---|
| `count_ctrl`, `n_chr_ctrl`, `freq_ctrl` | observed in the control group (`params.ControlGroup`) |
| `ci_lo_ctrl`, `ci_hi_ctrl` | 95 % Wilson interval — a normal approximation goes below zero for the rare alleles that are most of an HLA table |
| `count_case`, `n_chr_case`, `freq_case` | the cases, reported beside and **not** compared |
| `count_ref`, `n_chr_ref`, `freq_ref` | the reference. Null (N) alleles are dropped and the rest renormalised, because a null means the gene is absent from the haplotype — what this pipeline calls `not_typed` |
| `diff_ctrl_ref` | `freq_ctrl − freq_ref` |
| `fisher_p_ctrl_ref` | Fisher exact against the reference **counts**, unadjusted. The reference is an estimate too, and is treated as one — at 61,424 individuals its own interval is narrow, but the test still uses counts on both sides rather than a known *p* |
| `in_ref_only`, `in_obs_only` | present in one source and absent from the other |

`allele_frequency_summary.tsv`, one row per locus: `n_alleles`, both denominators,
`pearson_r` between the frequency vectors, `n_shared`, `max_abs_diff` and the allele
where it occurs, and the two directional counts.

**`n_shared` is the column that makes `r` readable**, and it is why there is no rank
correlation beside it. Most of each locus's vector is an allele one source carries and the
other does not, so a rank statistic over the whole vector is dominated by ties at zero and
its sign is set by the size of the two disjoint sets rather than by agreement. `n_shared`
says how many alleles the two sources actually have in common: 114 at HLA-A, but **2 at
HLA-E and 2 at HLA-F**, whose `r` of 0.999 and 1.000 is therefore computed on two points
and a mass of zeros. Read `r` and `n_shared` together or neither.

This compares **population frequencies**, not genotypes. It cannot say a given sample was
typed correctly. Of the 13 loci, 12 are comparable and 11 of those agree at *r* ≥ 0.98; DRB4
needs the null allele removed to get there. **DRB3 is excluded**, because the reference has no
way to record "this chromosome has no DRB3" and assigns one anyway — see
[OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) §3. The `comparable` and `note` columns of
`allele_frequency_summary.tsv` carry that decision and its reason. See
[METHODS.md](METHODS.md) §11.

## `allele_pgroup_map.tsv`

The crosswalk the association reads. One row per allele this cohort actually carries — all
2,298 of them, over all 33 typed genes, not only the 13 the reference covers:

| column | meaning |
|---|---|
| `gene`, `allele` | as they appear in `allele_dosage.tsv`'s column names |
| `p_group` | the IPD-IMGT P group, which is how the reference names the same allele |
| `in_reference` | 1 if the reference carries that P group, 0 if it has never observed it |
| `freq_reference` | its frequency in the reference, blank when `in_reference = 0` |

Two uses, both at reporting time and both one join away:

- an allele hit gets a **Japanese population frequency from 61,424 individuals** attached, in
  place of the 105 the previous reference could offer;
- **`in_reference = 0` is an artefact flag.** A hit on an allele that a 61,424-person Japanese
  panel has never observed is far more likely a typing error than a finding. 1,332 of the 2,298
  distinct alleles are flagged, but they are rare — chromosome-weighted the figure is 4.48 %, and
  that is the number to quote.

It is written by `ALLELE_FREQ_CHECK`, not `RESIDUE_MATRIX`, deliberately: see
[METHODS.md](METHODS.md) §13.

## `allele_frequency_sample.tsv`

One row per sample, the per-sample form of the same tail metric `typing_confound.tsv`
aggregates: `n_chr` typed chromosomes, `n_unconfirmed` of them carrying an allele absent
from the reference allele set, the resulting `frac_unconfirmed`, and `group`, `platform` and
`observed_depth` carried across from the manifest so the tail can be stratified any way.

A chromosome counts as unconfirmed when its call is `called` or `hemizygous` and its **P group**
is not in the reference set for that gene. It is **not** an error count: a real Japanese allele
that even a 61,424-person panel never sampled lands here too. It is only ever read as a
*relative* quantity, one stratum against another.

The cohort-wide figure fell from 10.8 % to **4.6 %** when the reference changed from 1000
Genomes JPT to jMorp on 2026-08-26. The typing did not change; the previous reference was small
enough that much of the apparent tail was simply an allele it had never seen.

`allele_frequency_check.tsv` **cannot** be used to derive this. It lists only alleles seen in
the CONTROLS or in the reference, so an allele seen only in cases is absent from it entirely
and reading the tail off that table silently undercounts the cases. This file, and
`typing_confound.tsv`, are computed against the reference allele set directly.

## `typing_confound.tsv`

The technical-artefact table, one row per stratum: the two phenotype groups, each platform,
and the depth-matched 17–21× window. Columns are `n`, `median_depth`, `n_chr`,
`n_unconfirmed`, `share_unconfirmed`, plus three tests filled in only where they apply —
`fisher_p_matched` on the matched window, `mannwhitney_p_depth` on cases against controls,
and `fisher_p_group` on the case/control contrast under each reference panel. The last two
rows are the SECONDARY panel's controls and cases, which is what makes the table self-
contained: the same contrast under two references, in one file.

The point of the matched window: the `15x` / `30x` in a platform label is a **target**, not a
measurement. Measured depth does not differ between cases and controls (18.65× vs 18.72×,
Mann–Whitney *P* = 0.194), so the case/control gap in the tail cannot be attributed to depth.
Restricted to samples measured at 17–21×, where the two medians are near 18.5×, the gap
persists (4.48 % vs 3.16 %, Fisher *P* = 3.58e-07) — it tracks the platform, not the
coverage.

The point of the second panel: "unconfirmed" is defined by a finite reference, so the level
could be an artefact of what that reference carries. Against the smaller panel the level
moves to 7.42 % against 5.40 % while the ratio moves only from 1.31 to 1.37. The level
belongs to the panel; the contrast does not.

Depth still matters *within* a platform: the tail runs 10.66 % at a median 14.17× down to
1.32 % at 30.93×. There is simply no case/control depth difference to act on. See
[OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) §3.

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
