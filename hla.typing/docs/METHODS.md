# METHODS — hla.typing

HLA allele typing and amino-acid residue calling. This document is about *why* each choice
was made; the output tree and column definitions are in [OUTPUTS.md](OUTPUTS.md), and this
study's own numbers are in [STUDY_NOTES.md](STUDY_NOTES.md).

---

## 1. What a call from this component means

An HLA-HD allele call is an estimate made from short reads against a dictionary of 46,005
known alleles. It is not a sequenced haplotype. Two things follow.

**Resolution is a result, not a constant.** A call may land at 2, 3 or 4 IMGT fields
depending on how much of the polymorphic sequence the reads separated. The depth achieved
is reported per sample and per platform, because a systematic difference in it between
groups would look exactly like a systematic difference in biology.

**A residue is derived, not observed.** The residue matrix is the allele call looked up in
IMGT's published protein alignment. Its correctness is entirely inherited from the allele
call; nothing here re-measures the protein.

## 2. The sample set is inherited, not decided here

`Flag_JHRPv6 == True` defines the v6 cohort — 3,592 samples. **Which of them is analysable is
decided upstream**, by `wgs.auto.par`'s sample QC, and this component reads that list rather
than repeating the judgement:

```
Flag_JHRPv6 == True                     3,592
  ∩ 07_sample_qc/…/sample_qc.keep.id    3,569      (23 excluded)
```

All 23 exclusions are heterozygosity outliers (`n_fail_het`), and that is the most relevant
exclusion this analysis could inherit: typing resolves a sample into **two haplotypes**, and a
sample whose genome-wide heterozygosity is anomalous is the one least likely to do so cleanly.
Typing them anyway would put the least trustworthy calls into the matrix with nothing marking
them.

If the keep list is absent the run **fails at launch naming the path it wanted**, rather than
silently typing 23 samples the rest of the study has already dropped.

## 3. Extraction: four sources of reads, and why none is optional

Reads are taken from

| | |
|---|---|
| `chr6:29,540,000–33,420,000` | the classical MHC |
| all `chr6_*_alt` scaffolds | divergent haplotypes that do not fit the primary assembly |
| all 525 `HLA*` contigs | hs38DH's HLA decoys — where reads from an allele far from the reference actually map |
| `*` | unmapped, where a read whose allele is in none of the above ends up |

Dropping any one of these loses reads from exactly the samples whose alleles are most
unusual, which is the opposite of what typing needs. The list is built **once per run** from
`hs38DH.fa.fai`; deriving it by parsing the 3 GB FASTA in every task costs about 200 s per
sample and buys nothing.

**The index is resolved by `realpath`.** 3,117 of this study's 3,569 CRAM paths are
symlinks whose `.crai` sits beside the real file, so the obvious `<cram>.crai` does not
exist and `samtools` fails with `[E::cram_index_load] Could not retrieve index file`. The
resolved path is passed with `-X`. A further 111 CRAMs have no index anywhere; those are
indexed **into this component's own tree**, never beside the source, which is a directory
we do not own.

**Reads are collated before conversion.** The extracted BAM is coordinate-sorted, so
`samtools fastq` without `collate` emits R1 and R2 in different orders *and different
counts*. HLA-HD survives that because it aligns the two files separately and re-pairs by
name — but the files are then wrong for every other consumer, silently. The pipeline
asserts R1 and R2 have equal counts.

## 4. Typing: HLA-HD, and the fact that it cannot report failure

`hlahd.sh -t 4 -m 100 -c 0.95` against the pinned dictionary and the pinned
`HLA_gene.split.3.50.0.txt`. That 1 KB gene list is pinned for the same reason the
dictionary is: it defines what every result is validated against, so a silent swap of it
in the group install would start failing samples whose typing was fine.

**HLA-HD always exits 0.** Its estimation step is a generated shell script with no `set -e`
that checks none of its children. When `hla_est` is killed — an out-of-memory kill, which
takes the most polymorphic genes first, so HLA-A and HLA-B — the tool still reports success
and a truncated result is published. In the reference implementation this happened to
**6 of 348 samples, reproducibly, with every exit code 0 and an empty stderr**.

The truncation is not visible downstream either. `pick_up_allele` writes
`Couldn't read result file.` on a line of its own, **carrying no gene name**, so a parser
reading by row position silently re-labels every gene after it.

Two mechanisms answer this, and **they only work together because the first runs inside
the second's process**:

- **`validate_typing.py`** checks that all 33 genes are present, that no unreadable-result
  token appears, that every `.est.txt` exists, and that no kill or OOM line is in the run
  log. Any one failing is a non-zero exit. It is the only signal there is.

  It is the **last line of the `HLAHD` process script**, not a process of its own. That
  matters: a downstream process cannot re-run an upstream one, so as a separate stage its
  failure could only terminate the run — and by then `publishDir` had already written the
  truncated result to `02.typing/`. Inside `HLAHD`, the same failure is a non-zero exit of
  the typing task, which is retried, and which publishes nothing.

- **A memory ladder**: 36,568 MB, retried at 73,136 MB. The reference used 2 GB, which is
  what produced the corruption. 36,568 is 4 × 9,142 — the whole memory four cores are
  already entitled to at this site — so the margin over the 3.6–5.3 GB measured on the
  pilot costs nothing.

`-t 4` rather than 1: `hlahd.sh` uses it for `bowtie2 -p` *and* for 4-way parallel
`hla_est`. Measured on one sample, that is 8 minutes against 30–50.

**One file goes downstream, not the directory it came from.** The validator looks at the
whole `result/` tree, and `HLAHD` publishes the whole tree, but the workflow passes on only
`<sample>_final.result.txt` — all `COLLECT_ALLELES` reads. Passing the directory would
stage 67 files per sample into that single downstream task: 239,123 symlinks for this
cohort against the 3,569 it needs.

Those files are staged under **HLA-HD's own names**, which are unique per sample, and
`collect_alleles.py` takes the sample id from the filename. The earlier design staged them
as `res_1 … res_N` and passed a parallel list of `SAMPLE=res_N` arguments, so two
independent orderings had to agree — and Nextflow does not number `res_*` at all when the
collection holds a single file, which made a one-sample run fail outright.

## 5. Concurrency is limited on I/O, not on compute

`EXTRACT_READS` and `INDEX_CRAM` are capped at `params.maxForksExtract` (48). Between them
they read ~88 TB of CRAM, and the filesystem they read it from is **shared** — running 3,569
of them as fast as the scheduler allows would degrade it for everyone, not only for this run.
Neither is CPU-hungry, so the cap costs little: the pilot took ~6 minutes per extraction, and
48 at a time keeps extraction off the critical path.

`HLAHD` is deliberately **not** capped. It is compute-bound, its memory is its real
constraint, and bounding compute is what the scheduler is for.

## 6. Residues: IMGT's alignment is read, never recomputed

IPD-IMGT/HLA publishes curated protein alignments. This component parses them. The
alternative — translate the dictionary's exon sequences, `blastx` them against a protein
database, and re-align with MUSCLE — costs three things at once, and all three were
measured in an implementation that does it:

1. **Positions stop being IMGT positions.** They become alignment column indexes that
   restart at 0 for every exon and shift with every release. A result at "DRB1 position 13"
   is then neither reproducible nor the position anyone else means by that name.
2. **`blastx` returning no hit yields no peptide**, and the obvious handling writes the
   string `"None"` into the reference as if it were a sequence. MUSCLE uppercases it and
   maps `o` to `X`, producing the four-residue "peptide" `N X N E`. In that implementation
   151 exon-groups are affected; `DQB1` exons 5 and 6, `DRB1` exon 5 and `DRB5` exons 5–6
   are **entirely fabricated**, and several of the affected groups are carried by every
   sample in the cohort.
3. **`blastx` SEG-masks low-complexity stretches**, replacing real residues with `X` —
   4.5 % of cells, concentrated in transmembrane regions.

None of those failure modes exists here, because no alignment is computed.

Parsing the alignment correctly needs three conventions, each verified against the file's
own `|` anchors:

- the alignment is written in **blocks**, each with its own `Prot` header. Only the first
  block's value is a starting point; taking the last one puts every position several
  hundred residues out;
- an allele whose leader is shorter than the alignment's is **padded with leading spaces**,
  and those are numbered positions where it has no residue. Stripping them shifts
  everything by the width of the padding — six positions for HLA-A, the difference between
  reading position 1 as `G` and reading it as `Y`;
- a `.` in the reference row is an insertion carried by another allele: it occupies a
  column but has no IMGT position. And IMGT numbering has no zero.

**Verified against known biology.** HLA-A positions 1–8 read `GSHSMRYF`, the canonical
mature N-terminus. *DRB1* at positions 11/13/71/74 gives `DRB1*04:05` = V/H/R/A,
`*15:01` = P/R/A/A, `*03:01` = S/S/K/R — the published values.

## 7. Genes: 33 typed, 19 with residues

HLA-HD types 33 genes. IPD-IMGT/HLA publishes a protein alignment for 19 of them:
`A B C DMA DMB DOA DOB DPA1 DPB1 DQA1 DQB1 DRA DRB1 DRB3 DRB4 DRB5 E F G`.

The other 14 — `DPA2 DRB2 DRB6 DRB7 DRB8 DRB9 H J K L T V W Y` — keep their allele calls
and are **recorded as having no residue reference**. An implementation that drops them
inside a bare `except` loses five genes and fourteen gene/exon pairs with no diagnostic.

## 8. Four outcomes, kept apart

A locus with no genotype is not one thing:

| | |
|---|---|
| `called` | two alleles |
| `hemizygous` | one allele and HLA-HD's `-` |
| `not_typed` | HLA-HD's `Not typed` — the locus is absent from this sample |
| `failed` | no row for the gene at all, or the unreadable-result token |

`failed` is a **pipeline** failure, not biology. Collapsing it into "missing" turns a
technical artefact into apparent absence, and if it is unbalanced between groups it is a
differential-missingness confounder with nothing in the record to say so.

## 9. Ambiguity is recorded, and resolved cohort-independently

When HLA-HD cannot choose it emits every candidate pair. Resolving those by their frequency
**in the current batch** makes a sample's genotype depend on who else was run with it: the
same sample genotypes differently in a pilot and in the full cohort, and the pipeline stops
being incrementally reproducible.

Here every candidate is written to `allele_ambiguity.tsv`, and the reported call is the
genotype at the deepest field depth every candidate **pair** agrees on. Pairs are compared
as *unordered* sets, because HLA-HD does not keep its two slots in a consistent order
between candidates — `G*01:04:03/G*01:01:08` and `G*01:01:01/G*01:04:01` are the same
genotype written twice. Comparing slot-wise collapses that to `G*01`; comparing as pairs
recovers `G*01:01 / G*01:04`.

## 10. Encoding: allele dosage, sorted diplotypes and 0/1/2 residue dosage

`residue_diplotype.tsv` sorts the two residues, so a heterozygote is one value. Writing
them in call order makes `LX` and `XL` different genotypes and splits every heterozygote
across two categories — in one measured implementation, 348 of 359 heterozygous codes had
both orientations present in the same file.

`residue_dosage.tsv` gives one 0/1/2 column per `(gene, position, residue)`. It is
order-free by construction and is what a regression takes.

Only positions polymorphic **in this cohort** are emitted. A monomorphic position carries
no information and cannot be tested.

`allele_dosage.tsv` is the same encoding one level up: one 0/1/2 column per **2-field
allele**, over every gene that produced calls, not only the 19 with residues. HLA
association is conventionally run at three levels — classical allele, amino-acid residue,
haplotype — and the allele test is the one reported first. It needs no reference at all,
being a re-encoding of `allele_calls.tsv`, but it is written by the same script as the
residue matrices so that the sample order is identical by construction rather than by
assumption.

**How each allele reached the reference is recorded once per allele**, in
`allele_match_map.tsv`: `allele`, `matched_to`, `how`, `n_chromosomes`. Recording it per
*(sample, allele)*, as an earlier version did, produces ~120,000 rows for this cohort of
which essentially all are ordinary successful prefix matches, and the few alleles that
match nothing are invisible among them. Those few, and only those, are in
`allele_unmatched.tsv`.

## 11. The external check, and the one it is not

Every other QC output here measures whether a call was **made** and how deeply it resolved.
None of them can see a call that was made confidently and is wrong.

`allele_freq_check.py` is the one check that can. The controls are an unselected Japanese
population sample, and the 1000 Genomes panel publishes 2-field types for **JPT**, the only
Japanese population in it — 105 samples over A, B, C, DQB1 and DRB1. Comparing the two
frequency spectra costs nothing and would move visibly if the typing were systematically
mis-assigning alleles.

Three things about it are load-bearing:

- **The reference is the 1KG panel, not HLA-HD's `freq_data/`.** The latter is the obvious
  candidate and is wrong twice: it is a *global* registry count (`A*01:01` = 326,922), not
  a Japanese one, and HLA-HD uses it to break its own ties.
- **The controls are compared, not the cases.** Cases are a disease series, and the MHC is
  the last place a disease series should be expected to match a population reference. Their
  frequencies are in the table beside, never as the comparison.
- **DQB1's reference denominator is 186 of 210 chromosomes**, in the reference itself. It
  is reported per locus rather than assumed to be 2*n*.

**This is not an accuracy measurement.** It compares population frequencies; it cannot say
a given sample was typed correctly, and a set of errors preserving the spectrum is
invisible to it. Only typing samples with known types measures accuracy. That was scoped —
104 JPT samples have both a published type and a downloadable 30× CRAM — and deferred; see
[OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) §2 for the cost and the decision. Truncating 150 bp
samples to 100 bp and re-typing the same individuals, which would give the read-length cost
directly and with no confounding, is recorded there too.
