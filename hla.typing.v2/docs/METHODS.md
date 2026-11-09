# hla.typing.v2 — methods

**Reported evidence.** The basic claim that the genotypes are usable rests on two checks: completeness
(§2–3) and the agreement of control allele frequencies with ToMMo 61KJPN (§5, r). Everything else here —
read evidence, confidence classes, *DRB* reconciliation, 1KG JPT, homozygosity, tiers — is kept as
supplementary material.

## 1. Scope

This component has two aims:

1. **HLA genotypes for the full_mainland cohort** (3,101 samples: 2,662 AGP3K controls and 439 PH cases).
2. **Evidence of their quality.**

**The typing itself is not re-run.** HLA-HD ran once over these 3,101 samples in `hla.typing` (2026-08-25).
This component starts from that output and derives everything else from it.

**Out of scope.**
- **Platform and sequencing-depth differences.** These are left to a later down-sampling and
  platform-sensitivity study; `platform` and `observed_depth` are carried as per-sample columns only.
- **Association.** The matrices the association reads are produced (§7 exports, §8 association inputs, one set per cohort); no test is run here.

**What the quality evidence can and cannot show.** There is no per-sample truth here: no samples with
laboratory HLA types were typed. The evidence can therefore establish that the calls are consistent with
the reads behind them and with the Japanese population. It cannot establish per-sample accuracy, and it is
blind to errors that leave allele frequencies unchanged.

## 2. Inputs: the typing run, archived and checked

**Reads.**
- Regions: chr6:29,540,000–33,420,000, the 16 chr6 alternate contigs, the 525 HLA contigs of hs38DH and the
  unmapped reads — 543 regions.
- They were extracted from each CRAM with `samtools view`, then `collate`, then `fastq`.
- Secondary and supplementary records were dropped; duplicates and QC-fail reads were kept.

**Typing.** HLA-HD 1.7.1 (Kawaguchi et al., 2017):
- `hlahd.sh -t 4 -m 100 -c 0.95 -f freq_data`;
- dictionary IMGT-HLA-3.64.0_N150; gene split file 3.50.0, which lists 33 genes.

HLA-HD aligns the reads to the exon and intron sequences of every named allele and scores every allele pair
by the reads each allele explains, exon by exon. It reports the best pair. Where pairs tie, or an allele has
several candidates, it uses its built-in allele-frequency table to choose; only *A, B, C, DPA1, DPB1, DQA1,
DQB1* and *DRB1* have one.

**Archive (step 00a).** For every sample one `tar.gz` holds:
- `final.result.txt`;
- the 33 `<gene>.est.txt` (candidate pairs) and 33 `<gene>.read.txt` (reads per allele);
- the HLA-HD run log and its validation record;
- the HLA-HD and extraction task scripts and logs.

That is 75 files per sample, 167 MB in all. The originals are symlinks into `hla.typing/work/`, which
also holds about 0.8 GB per sample of intermediate alignments that are not needed here. The archive makes
this component independent of that tree.

**Provenance (step 00b).** The run is accepted as one consistent typing only if all of the following hold.
1. Every sample's run log prints the same HLA-HD version (1.7.1).
2. Every sample ran the same HLA-HD command and the same extraction command (sample name masked).
3. All 3,101 validation records read `ok`. That record means: all 33 genes present, every estimate file
   present, and no out-of-memory kill.

The extraction log gives the read pairs given to HLA-HD: median 431,582, 1st percentile 234,238. HLA-HD's
frequency table is copied into `results/00.inputs/hlahd_freq_data/`; it is unchanged since the run.

## 3. From HLA-HD output to genotypes (step 01)

**3.1 Candidates.** `<gene>.est.txt` lists the best-scoring pairs (`#Best allele pair` = the number that
tie). Each side of a pair is a list of allele names the reads cannot tell apart. The result line in
`final.result.txt` only displays this: when a side's candidates differ at the 2nd or 3rd field it prints
every combination. It is therefore used only to identify the pair HLA-HD reported, never to read alleles.

**3.2 Agreed name.**
- **Pairs are aligned first:** each tied pair is oriented to the reported pair, by the larger shared prefix,
  because HLA-HD does not keep the two sides in a fixed order.
- **Agreed name of each allele:** the deepest field on which all of its candidates agree, across all tied
  pairs.
- **Expression suffixes** (N, L, S, Q) are kept only when every candidate carries the same one.

**3.3 The 2-field allele.** Its source is the first of these rules that applies.

| source | rule |
|---|---|
| agreed | the candidates agree to at least 2 fields |
| prior | they do not, and HLA-HD printed one allele for them after applying its frequency table (8 genes). The ratio between the table's count for the printed 2-field allele and for the best other candidate is kept. |
| P group | the candidates all lie in one IPD-IMGT P group (identical protein over exons 2–3 for class I, exon 2 for class II; Anthony Nolan Research Institute, n.d.-b); named by the group's 2-field representative |
| none | no 2-field allele is written — never an arbitrary pick |

**3.4 States.** For the obligate genes (all reported genes except *DRB3/4/5*):

| state | when |
|---|---|
| heterozygous | two different 2-field alleles |
| homozygous | one allele (`-` in HLA-HD's output) or the same 2-field allele twice |
| ambiguous | agreement at 1 field only |
| unresolved | no field agreed |
| failed | no candidate, or no result line |

**3.5 *DRB3, DRB4, DRB5* copy number.** The number of *DRB* genes differs between DR haplotypes (Andersson,
1998), and the *DRB1* group determines which other *DRB* genes a haplotype carries, with occasional exceptions
(Anthony Nolan Research Institute, n.d.-a). The expected copies come from the two *DRB1* alleles' first
field:

| *DRB1* | haplotype group | paralogue |
|---|---|---|
| \*15, \*16 | DR51 | *DRB5* |
| \*03, \*11, \*12, \*13, \*14 | DR52 | *DRB3* |
| \*04, \*07, \*09 | DR53 | *DRB4* |
| \*01, \*08, \*10 | DR1, DR8, DR10 | none |

How each combination of expected and reported copies is resolved:

| expected | HLA-HD reports | result |
|---|---|---|
| 0 | nothing | absent |
| 0 | an allele | absent, flagged `drb_unexpected` |
| 1 | one allele | present, one copy |
| 1 | two alleles | the allele with more reads unique to it is kept, flagged `drb_reconciled` |
| 2 | one allele | present, two copies of it |
| ≥ 1 | nothing | missing (no call), flagged `drb_missing` |

- **HLA-HD's own call is always kept** beside the reconciled one.
- **When *DRB1* is unresolved**, HLA-HD's call is kept and flagged `drb1_unknown`.

## 4. Read evidence and call confidence (step 02)

**Reads.** From `<gene>.read.txt`:
- the reads assigned to each allele;
- the reads assigned to one allele but not the other — **unique reads**.

Shared reads fit both alleles and say nothing about whether a heterozygote really carries two alleles. The
unique reads cover the positions where the two alleles differ, and they are that evidence.

**Flags.** Thresholds were fixed in `config.py` before this step was run.

| flag | rule |
|---|---|
| `low_support` | an allele with fewer than 10 reads |
| `imbalanced` | a heterozygote whose minor allele has less than 15 % of the unique reads; assessed only with at least 10 unique reads, because closely related alleles leave few |
| `weak_prior` | a prior-chosen allele whose count ratio is below 100 |
| `pgroup` | a P-group resolution |
| `drb_*` | from §3.5 |

**Confidence.**

| class | rule |
|---|---|
| high | 2-field alleles agreed by the reads, or chosen by the prior against alternatives at least 100 times rarer, and no flag |
| medium | exactly one weakness |
| low | two or more weaknesses, or no 2-field allele (ambiguous) |
| no call | failed, unresolved, missing |

**Sample QC.** A sample is listed as an outlier if any of these holds:
- read pairs below the cohort's 1st percentile;
- any failed obligate gene;
- three or more reported genes without a 2-field allele.

Outliers are not removed.

## 5. External panels (step 03)

**Panels.**
- **ToMMo 61KJPN:** allele frequencies of 61,424 Japanese individuals at 13 loci (Tadaka et al., 2024),
  named by IPD-IMGT P group.
- **1000 Genomes JPT:** 105 Japanese individuals at 5 loci, typed in silico from the 1000 Genomes data
  (Abi-Rached et al., 2018).

**Common rules.**
- **Only controls are compared.** Cases are a disease-selected series whose HLA frequencies may differ for
  real reasons; *HLA-B\*52:01* and *DPB1\*02:02* have been reported to be more frequent in Japanese CTEPH
  (Tanabe et al., 2005; Kominami et al., 2009).
- **All comparisons are made on P groups.**
- **Null alleles are dropped** on both sides and the rest renormalised. ToMMo writes *DRB4\*03:01N* for most
  chromosomes without *DRB4*, but not cleanly: its non-null *DRB4* share (45.7 %) exceeds the DR53 share of its
  own *DRB1* (35.3 %).
- **Who contributes:** a control contributes to a locus only if every copy it carries has a 2-field
  allele.

**Loci whose frequencies are compared: 11** (*A, B, C, E, F, G, DRB1, DQA1, DQB1, DPA1, DPB1*). Pearson r is
computed within each locus, never across loci. The other eight reported genes are not compared, each for a
stated reason (the rows of *DRB3* and *DRB4* stay in `locus_concordance.tsv` with `frequency_comparable = 0`
and the reason; both still enter the absent-allele check):
- ***DRA, DMA, DMB, DOA, DOB, DRB5*: no panel.** ToMMo publishes no frequencies for them, and 1KG JPT covers
  *A, B, C, DRB1, DQB1* only.
- ***DRB3*: not comparable.** ToMMo assigns a *DRB3* allele to every chromosome (*DRB3\*01:01P* = 0.674),
  although only DR52 haplotypes carry the gene (25.4 % of chromosomes by ToMMo's own *DRB1*).
- ***DRB4*: not compared (decided 2026-10-20).** Its presence follows *DRB1*, which is compared already;
  ToMMo's *DRB4* entries do not match the DR53 share of its own *DRB1* (45.7 % vs 35.3 %); and among carriers
  most alleles are *DRB4\*01:01* (92 % here, 99 % in ToMMo), so r carries little information. Until then it was
  counted, with r = 1.000.

The DR haplotype groups themselves agree: by *DRB1*, DR51 / DR52 / DR53 make up 20.7 / 25.4 / 35.3 % of ToMMo
chromosomes and 19.0 / 24.6 / 33.5 % of our controls' (from the copies of *DRB5 / DRB3 / DRB4*). This is
*DRB1* evidence, not evidence for *DRB3/4/5*, because our copies are set from *DRB1* (§3.5).

**Two features of ToMMo that the metrics must allow for.**
1. ***DRB3* frequencies are not comparable** (above).
2. **Dictionary vintage.** The newest allele in ToMMo was first named in IPD-IMGT/HLA release 3.21.0; HLA-HD
   here used 3.64.0. An allele first named after the panel's newest allele cannot appear in the panel
   whatever its true frequency. Each allele is dated with the IPD-IMGT allele history file (release 3.64.0,
   pinned in `assets/imgt/`).

**Metrics, per panel and locus.**

| metric | definition |
|---|---|
| r | Pearson correlation of allele frequencies |
| D | ½ Σ \|f<sub>ours</sub> − f<sub>panel</sub>\|: the smallest share of chromosomes that would have to change allele to make the two distributions equal |
| D<sub>null</sub> | D expected from sampling alone: both sides redrawn 1,000 times at their own chromosome counts from their pooled frequencies |
| excess D | D − mean D<sub>null</sub> |
| D<sub>shared</sub> | D over the alleles both sides carry, each renormalised to them; also computed on high-confidence calls only |
| absent (assessable) | share of our chromosomes on an allele the panel does not carry although its dictionary could contain it |
| absent (after vintage) | the same for alleles named after the panel's newest allele; reported, not judged |
| homozygosity | observed vs Hardy–Weinberg expectation from our own 2-field allele counts (Levene's exact expectation); 95 % Wilson CI; permutation P (1,000 shuffles) |
| dropout | (observed − expected) / (1 − expected): the share of heterozygotes reported as homozygotes, if Hardy–Weinberg holds; with the CI from the observed proportion |
| shift direction | for each allele whose frequency differs from ToMMo by at least 1 point at the five 1KG loci: whether 1KG JPT lies on our side or on ToMMo's |
| panel ceiling | D between ToMMo and 1KG JPT, beside ours against 1KG JPT |
| confidence check | the absent-allele share by call confidence (§4) |

## 6. Locus tiers (step 04)

**Criteria.**

| criterion | met when | applies to |
|---|---|---|
| absent alleles | assessable share ≤ 2 % | ToMMo loci |
| dropout | lower 95 % bound ≤ 3 % | obligate genes |
| call quality | high-confidence calls ≥ 95 % | every reported gene |

**Rule.** Tier = 1 + the number of criteria failed, capped at 3. A locus with no external panel (*DRA, DMA,
DMB, DOA, DOB, DRB5*) starts at 2, because nothing outside our own data speaks for it.

**Revision record (2026-10-19).** The rule written before step 03 first ran had four criteria:
- excess D against ToMMo ≤ 0.03;
- absent alleles ≤ 2 %, counting every absent allele;
- the CI of observed/expected homozygosity including 1;
- high-confidence calls ≥ 95 %.

Three changes were made after the first run, each because the criterion did not measure what it was meant
to.

1. **Frequency discordance moved from tiering to description.**
   - Excess D against ToMMo was 0.04–0.11 at *A, B, C, DRB1, DQB1, DPB1, G* and *DQA1*.
   - Its part over the alleles both sources carry did not shrink when only high-confidence calls were used
     (*A* 0.085 → 0.078, *B* 0.112 → 0.121).
   - For 24 of the 33 alleles that differ from ToMMo by at least 1 point, the independent 1KG JPT sample lies
     on our side.

   The shifts are therefore consistent with differences between cohorts rather than with typing errors,
   and cannot serve as an error rate.
2. **Alleles named after ToMMo's dictionary are no longer counted** in the absent-allele share: 880 of the
   1,355 absent chromosomes (65 %). They are reported as a separate column.
3. **The homozygosity criterion became the dropout criterion with a 3 % tolerance.** The original test
   tolerated no excess at all; with 2,600 controls, almost every locus fails it. 3 % is the tolerance
   already used for the frequency criterion.

The thresholds that were kept (2 %, 95 %) were not changed.

## 7. Exports (step 05)

These keep the file names and layouts of `hla.typing`:
- `allele_dosage.tsv`: sample × 2-field allele, 0/1/2 copies;
- `residue_dosage.tsv`: sample × `gene:IMGT position:residue`, 0/1/2 copies;
- `allele_pgroup_map.tsv`;
- plus `call_confidence.tsv`, so that calls can be filtered.

What differs from `hla.typing`'s exports:
- only the 19 reported genes;
- *DRB3/4/5* follow §3.5, so a single copy is dosage 1 and an absent gene is 0;
- a call without a 2-field allele is blank;
- **a residue the alignment does not determine is blank.** Residues come from IPD-IMGT's own protein
  alignments (release 3.64.0), read by `build_residue_ref.py` (copied unchanged from `hla.typing`); each
  2-field allele is matched to the alignment exactly, else to its first deeper allele (same two fields, so
  the same protein); never to a shorter name, which would borrow another protein. One called allele has
  neither — *B\*44:466N*, a null allele with one copy, absent from the 3.64.0 alignment — so its carrier is
  blank at every *B* residue (`allele_unmatched.tsv`). Where one of a
  sample's alleles has no determined residue at a position (`*` unsequenced, `.` gap), the sample is blank
  at every residue of that position. Written as 0, it would claim the sample carries none of them. A
  position is exported when at least two residues are determined among the samples, one
  column per residue seen; so at every called position of an obligate gene the dosages sum to 2.

## 8. Association inputs (step 06)

One set of matrices per cohort of the association (`config.COHORTS`: narrow, intermediate and full
mainland, from the phenotype files of `wgs.auto.par`). Every sample of a phenotype file must be typed;
none is dropped. Each exported column is measured within the cohort and kept or dropped by these rules,
in this order; `markers.tsv` records the first one it fails:

| rule | a column is dropped when | why |
|---|---|---|
| varies in the cohort | it has no copy of its less common form | it cannot inform a test |
| genes of the association | it belongs to *DRB3*, *DRB4* or *DRB5* | presence follows *DRB1* (§3.5); no GRCh38 primary locus for *DRB3/4*; ToMMo cannot check them |
| call rate ≥ 95 % in cases and in controls | under 95 % of the cases, or of the controls, have a genotype at the marker | completeness, measured in each group because sequencing platform and depth differ between them; set below every gene's completeness in either group (lowest *G* in controls, 97.8 %) |
| MAF ≥ 0.5 % | the less common form is under 0.5 % of the copies typed | the floor of comparable HLA studies; ~4 copies expected in cases at the floor |

Result (narrow / intermediate / full): 143 / 144 / 142 allele and 672 / 681 / 679 residue columns.
Call rate removes no allele column, because every allele of a gene shares that gene's completeness.
Measured overall instead, the call-rate rule would keep the same columns in all three cohorts (lowest group
call rate among the tested columns: 96.7 % in cases, 95.4 % in controls); the per-group rule guards later
re-runs. `verify.sh` recomputes both groups' call rates and MAF from the written matrices and checks their md5. **No test is
run here; the association analysis reads these matrices as they stand.**

## 9. Figures (step 07)

Eight figures, journal width, catalogued in `results/figures/README.md`. They use the project's shared
style (`analysis/_shared/scripts/plot_style.py`, imported unchanged) and Liberation Sans (Arial metrics,
bundled in `assets/fonts/`).

## 10. Limitations

- **No per-sample accuracy.** Errors that preserve frequencies are invisible.
- **Read counts include PCR and optical duplicates,** which were not removed at extraction, so support and
  balance count reads, not molecules.
- **Most absent alleles are not caught by read-level evidence.** HLA-HD sometimes chooses a pair of rare
  alleles that explains the reads as well as a common pair, with balanced unique reads (e.g. *B\*40:06 +
  B\*40:74*). Such calls are high-confidence by read evidence; 293 of the 475 assessable absent chromosomes
  are of this kind.
- **ToMMo's typing method is not stated in the pinned copy of the panel.** If it shares HLA-HD's biases,
  frequency agreement cannot reveal them.

## References

Abi-Rached, L., Gouret, P., Yeh, J. H., Di Cristofaro, J., Pontarotti, P., Picard, C., & Paganini, J. (2018).
Immune diversity sheds light on missing variation in worldwide genetic diversity panels. *PLoS ONE, 13*(10),
e0206512.

Andersson, G. (1998). Evolution of the human HLA-DR region. *Frontiers in Bioscience, 3*, d739–d745.

Anthony Nolan Research Institute. (n.d.-a). *DRB haplotypes*. HLA Nomenclature. Retrieved October 19, 2026,
from https://hla.alleles.org/pages/other/drb_haplotypes/

Anthony Nolan Research Institute. (n.d.-b). *P groups*. HLA Nomenclature. Retrieved October 19, 2026, from
https://hla.alleles.org/pages/p&g_groups/p_groups/

Barker, D. J., Maccari, G., Georgiou, X., Cooper, M. A., Flicek, P., Robinson, J., & Marsh, S. G. E. (2023).
The IPD-IMGT/HLA Database. *Nucleic Acids Research, 51*(D1), D1053–D1060.

Kawaguchi, S., Higasa, K., Shimizu, M., Yamada, R., & Matsuda, F. (2017). HLA-HD: An accurate HLA typing
algorithm for next-generation sequencing data. *Human Mutation, 38*(7), 788–797.

Kominami, S., Tanabe, N., Ota, M., Naruse, T. K., Katsuyama, Y., Nakanishi, N., ... & Kimura, A. (2009).
HLA-DPB1 and NFKBIL1 may confer the susceptibility to chronic thromboembolic pulmonary hypertension in the
absence of deep vein thrombosis. *Journal of Human Genetics, 54*(2), 108–114.

Tadaka, S., Kawashima, J., Hishinuma, E., Saito, S., Okamura, Y., Otsuki, A., ... & Kinoshita, K. (2024).
jMorp: Japanese Multi-Omics Reference Panel update report 2023. *Nucleic Acids Research, 52*(D1), D622–D632.

Tanabe, N., Kimura, A., Amano, S., Okada, O., Kasahara, Y., Tatsumi, K., ... & Kuriyama, T. (2005).
Association of clinical features with HLA in chronic pulmonary thromboembolism. *European Respiratory
Journal, 25*(1), 131–138.
