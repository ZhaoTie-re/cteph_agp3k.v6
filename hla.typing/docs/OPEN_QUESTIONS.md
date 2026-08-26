# OPEN QUESTIONS — hla.typing

Known defects and unvalidated assumptions, deliberately left in place. Each entry says what is
wrong, what the evidence is, what the current behaviour is, and what would fix it. Nothing here
is a plan; these are decisions that have been *deferred*, and the deferral is the record.

Method and rationale for what the component *does* do: [METHODS.md](METHODS.md).

---

## 1. DRB3/4/5 hemizygotes are encoded as homozygotes

**Status: known-wrong, deliberately unchanged.** The current behaviour matches the reference
implementation this component was built alongside
(`/LARGE0/gr10478/home/yuxun.zhang/project_typing`), and changing one without the other would
make the two incomparable. Fix both, or neither.

### What happens now

HLA-HD writes `-` for a second allele it did not call. `collect_alleles.py` fills it with the
first:

```python
elif len(a) >= 2 and a[1] == '-':
    state, a1, a2 = 'hemizygous', a[0], a[0]
```

The reference does exactly the same, in `do0.HLA_allele_summary.py` (twice, lines 42 and 58):

```python
if str(r.iloc[l, 2*k+2]) == "-":
    r.iloc[l, 2*k+2] = r.iloc[l, 2*k+1]
```

Neither treats DRB3/4/5 differently from any other gene.

### Why that is wrong for these three genes and right for the others

`-` carries two different meanings, and nothing in the call distinguishes them:

- **At the classical loci it means homozygous.** HLA-A shows `-` on 13.5 % of samples, and
  13.5 % of people do not have one copy of HLA-A. Filling it is correct.
- **At DRB3/4/5 it means hemizygous.** These are haplotype-dependent paralogues that are
  genuinely absent from many chromosomes: DRB3 travels with `DRB1*03/*11/*12/*13/*14`, DRB4
  with `*04/*07/*09`, DRB5 with `*15/*16`, and `DRB1*01/*08/*10` carry none.

Copy number is therefore **determined by the DRB1 genotype**, which this component already
calls. Measured over the reference cohort's 232 samples:

| | DRB1 implies 0 copies | implies 1 | implies 2 |
|---|---|---|---|
| **DRB3** — HLA-HD says `Not typed` | 106 | 0 | 0 |
| HLA-HD says one allele + `-` | 6 | **99** | 5 |
| HLA-HD says two alleles | 0 | 6 | 10 |
| **DRB4** — `Not typed` | 67 | 0 | 0 |
| one allele + `-` | 0 | **117** | 21 |
| two alleles | 0 | 14 | 13 |
| **DRB5** — `Not typed` | 163 | 0 | 0 |
| one allele + `-` | 0 | **61** | 3 |
| two alleles | 0 | 3 | 2 |

`Not typed` corresponds to 0 copies with no exceptions. **One allele + `-` corresponds to 1
copy** — 99/110, 117/138 and 61/64 of the time.

### The consequence

A sample with one DRB4 copy gets dosage **2** for every DRB4 residue instead of 1. On the
reference cohort that is roughly 280 of 348 samples across the three genes.

**So the DRB3, DRB4 and DRB5 columns of `residue_dosage.tsv` AND of `allele_dosage.tsv` are
not trustworthy as dosages.** The defect is in `collect_alleles.py`, which writes the single
allele into both slots, so everything derived from `allele_calls.tsv` inherits it. The allele
calls themselves are fine; it is the copy number attached to them that is wrong. The other 16
genes are unaffected, and the frequency check of §2 does not touch DRB3/4/5 at all.

### What a fix would look like

Infer copy number from the DRB1 call, which needs no new data:

| inferred copies | encoding |
|---|---|
| 0 | dosage 0 at every residue — the gene is absent, which is not the same as a missing call |
| 1 | one allele, dosage 1 |
| 2 | the current diploid handling |

The change belongs in `collect_alleles.py`, at the `a[1] == '-'` branch that currently writes
`a1` into both slots — fixing it there fixes `allele_dosage.tsv` and `residue_dosage.tsv`
together, since both read `allele_calls.tsv`. The disagreeing minority (6 DRB3 samples called despite DRB1
implying 0 copies) should be recorded, not silently overridden either way.

The haplotype linkage is textbook, so this is not a novel method — but it does make this
component's DRB3/4/5 output differ from the reference's, which is the reason it is deferred.

---

## 2. Typing accuracy has never been measured per sample

Every QC number this component produces is about whether a call was *made* and how deeply it
resolved. **Nothing measures whether a given sample's call is right.**

### What was added instead, and what it does establish

`allele_freq_check.py` compares the controls' 2-field allele frequencies at A, B, C, DQB1
and DRB1 against the **1000 Genomes JPT** panel — the only Japanese population in it.
That is a real external check and it costs nothing: if the typing were systematically
mis-assigning alleles the frequency spectrum would move, and no call-rate metric would
notice. `figures/allele_frequency.png` draws it.

It is **not** an accuracy measurement. It compares population frequencies; it cannot say a
given sample was typed correctly, and a set of errors that happens to preserve the spectrum
is invisible to it. The reference is 105 samples, so the standard error on a frequency near
0.4 is ~0.034 — it bounds gross error and nothing finer. And it covers five loci: the
panel carries no DPB1, no DQA1 and no DRB3/4/5, so **the loci with the least reliable
typing here are exactly the ones it cannot check**.

### The measurement that would settle it, scoped and deferred

Typing samples whose types are already published. Scoped on 2026-08-25 and **deferred by
decision, not by cost**:

| | |
|---|---|
| truth set | `1KG_HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt` — **pinned in this project**, 2,693 samples, 2-field |
| usable samples | **JPT 105** have a published type; **104** of those have a downloadable 30× CRAM |
| why JPT and not the local CRAMs | the 80 1KG CRAMs already on this system contain **zero** JPT. `A*24:02` (0.39 here), `DRB1*04:05` and `B*52:01` are common in Japan and rare in CEU/YRI, so a non-Japanese truth set answers a different question |
| download | ~15.1 GB per CRAM → **~1.57 TB**. `1000genomes.s3.amazonaws.com` served 10.4 MB/s single-stream; EBI's FTP mirror only 115 KB/s. The `.crai` sits beside each CRAM |
| where the URLs are | `1KG_HLA_types/1000G_2504_high_coverage.sequence.index`, filtered on `POPULATION == JPT` |

Everything except the CRAMs is already on disk, so re-opening this is a download and a
`--MaxSamples`-style run, not a re-investigation.

The second experiment is unchanged and needs no download at all: **hard-trim 150 bp samples
to 100 bp and re-type the same individuals.** One variable, no confounding. The pilot's
per-platform comparison bounds this effect but cannot isolate it.

Until the first is done, per-sample accuracy rests on internal consistency plus a
population-frequency check.

---

## 3. 13-18 % of A/B/C chromosomes are called as alleles JPT never carries

**Found by the frequency check of §2 on the full 3,569-sample run — this is the defect that
check was built to detect, and it is the largest open question here.**

Per locus, the share of control chromosomes assigned to a 2-field allele that does not occur
once in the 1000 Genomes JPT panel, beside the total shortfall across alleles the panel has
at frequency >= 0.05:

| locus | tail alleles | tail chromosomes | share | common-allele shortfall |
|---|---|---|---|---|
| A | 413 | 861 | **13.8 %** | 0.131 |
| B | 418 | 1,107 | **17.8 %** | 0.143 |
| C | 347 | 797 | **12.8 %** | 0.126 |
| DRB1 | 152 | 434 | 7.0 % | 0.049 |
| DQB1 | 117 | 155 | 2.5 % | 0.019 |

**The two columns match.** The tail is not extra diversity sitting beside a correct common
spectrum; it is taken *out of* the common alleles, which is why 19 of the 25 alleles the
panel carries at >= 0.08 come out low in our controls — `A*24:02` 0.309 against 0.391,
`B*40:02` 0.063 against 0.101, `C*07:02` 0.108 against 0.152.

### What it is not

- **Not the ambiguity resolution of §5.** Traced per allele: `B*44:302` appears 67 times and
  exactly 1 of those came from an ambiguous call; `B*40:379`, `B*07:381`, `B*35:380`,
  `DRB1*14:253` and `A*24:50` came from ambiguity **zero** times. HLA-HD reports these as
  unambiguous calls.
- **Not the reference's small size.** 210 chromosomes cannot sample an allele below ~0.5 %,
  so a long tail of *names* is expected. A tail carrying 13-18 % of the **mass** is not:
  a Japanese population puts >95 % of chromosomes at these loci on ~30 common alleles.
- **Not a parsing artefact.** The reference's own `*`-flagged and `/`-ambiguous values were a
  separate bug, fixed; it moved `n_ref_only` from 11 to 1 and left these numbers intact.

### What it probably is

HLA-HD choosing a nearly identical rare allele over the common one out of a 46,005-allele
dictionary — `B*44:302:01` instead of `B*44:03`, `B*40:379` instead of `B*40:02`. The
frequency table meant to break exactly these ties (`-f freq_data`) **is** passed, and it is a
global registry count rather than a Japanese one, which is the wrong prior for this cohort.

### Why it is not fixed here

Every candidate fix changes genotypes, and the choice is a study decision, not a code one:

1. restrict the dictionary to alleles observed in East Asians, and re-type;
2. post-hoc collapse tail alleles onto their nearest common relative — cheap, reversible, and
   an invention nobody else uses;
3. accept it, and treat the 2-field allele test as the primary analysis only for alleles the
   reference panel confirms.

### IT IS DIFFERENTIAL BY PHENOTYPE, and that is the part that matters

Call rate is flat between the groups (0.917 vs 0.916), which is what made this look benign at
first. **Typing accuracy is not.** The share of chromosomes on a JPT-confirmed allele:

| locus | controls | cases | difference | Fisher *P* |
|---|---|---|---|---|
| A | 86.2 % | 90.6 % | **+4.4 %** | 2 × 10⁻⁴ |
| B | 82.2 % | 86.3 % | **+4.0 %** | 0.0025 |
| DRB1 | 93.0 % | 95.5 % | +2.4 % | 0.0053 |
| C | 87.2 % | 89.4 % | +2.2 % | 0.067 |
| DQB1 | 97.5 % | 98.5 % | +0.9 % | 0.10 |
| **all five** | **89.2 %** | **92.0 %** | **+2.8 %** | |

The cause is depth, and depth is the study's confound: **every control is HiSeqX 15×**, most
cases are 30×. Tail share by platform — NovaSeq 30× 2.3 %, DNBSeq-G400RS 30× 4.3 %,
DNBSeq-T7 30× 6.3 %, **HiSeqX 15× 14.8 %, DNBseq-G400RS 15× 20.4 %**; against depth,
Spearman ρ = −0.109, *P* = 8 × 10⁻¹¹.

So the controls are typed measurably worse than the cases, for a technical reason that cannot
be separated from phenotype. A common allele depleted more in controls than in cases reads as
enrichment in cases — a false risk association, from the typing alone.

*(An earlier version of this section said the defect was not differential. That was read off
call rate, which does not measure accuracy, and it was wrong.)*

### Which output survives it, measured

Comparing each tail allele's residues against the commonest confirmed allele of its 1-field
family:

| locus | pairs | median residues differing | share differing by ≤ 2 | positions |
|---|---|---|---|---|
| A | 396 | **2** | 51 % | 390 |
| C | 338 | **2** | 55 % | 417 |
| B | 395 | **5** | 41 % | 483 |
| DQB1 | 117 | 47 | 15 % | 274 |
| DRB1 | 144 | 85 | 22 % | 290 |

**Class I residues are largely insulated**: a mis-call to a rare allele of the same family
moves 2–5 of ~400 positions. **Class II residues are not** — DQB1 and DRB1 tail calls move
tens of positions — but their tails are also the smallest (2.5 % and 7.0 %).

Order of robustness, most to least: class I residues → class II residues → `allele_dosage.tsv`.
Rare rows of `allele_dosage.tsv` should not be tested without checking them against
`05.qc/allele_frequency_check.tsv` first.

---

## 4. The IMGT alignment parser is validated on three anchors

`build_residue_ref.py` parses IPD-IMGT/HLA's `*_prot.txt` files from scratch — block starts,
leading-space padding, `.` insertion columns, the no-zero rule. It is checked against:

- HLA-A positions 1–8 = `GSHSMRYF`, the canonical mature N-terminus
- *DRB1* `*04:05`/`*15:01`/`*03:01` at positions 11/13/71/74, matching the published values
- zero occurrences of the `NXNE` signature that a failed blastx round-trip produces

Three anchors on two genes is not a test suite. Note in particular that **DRB3, DRB4 and DRB5
all use `DRB1*01:01:01:01` as their alignment reference** — a cross-gene reference, which the
parser handles by construction but which was never verified against a known DRB3/4/5 residue.

---

## 5. Ambiguity resolution is a local invention

When HLA-HD cannot choose between candidate pairs, `collect_alleles.py` reports the genotype at
the deepest field depth at which every candidate *pair* agrees, comparing pairs as unordered
sets. It is deterministic and cohort-independent, which the reference's frequency-based
resolution is not — but it is not a published method, and no one else does it this way.
`allele_ambiguity.tsv` keeps every candidate, so the decision is reversible.

---

## 6. `residue_dosage.tsv` bakes in a filter

Only positions polymorphic **in this cohort** are emitted. That is a reasonable default and it
is stated, but it puts an analysis decision inside data production: an analyst cannot recover
the monomorphic positions without re-running the stage.
