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

`allele_freq_check.py` compares the controls' allele frequencies against **jMorp
61KJPN-HLA** — ToMMo's panel of **61,424 Japanese individuals over 13 loci**. That is a real
external check and it costs nothing: if the typing were systematically mis-assigning alleles
the frequency spectrum would move, and no call-rate metric would notice.
`figures/allele_frequency.png` draws it.

**Measured, 2026-08-26.** Eleven of the thirteen loci agree at Pearson *r* = 0.986–1.000:

| | | | | | |
|---|---|---|---|---|---|
| A | 0.996 | DPA1 | 1.000 | DQB1 | 0.992 |
| B | 0.986 | DPB1 | 0.998 | DRB1 | 0.988 |
| C | 0.992 | DQA1 | 0.997 | E · F · G | 0.999 · 1.000 · 0.993 |
| **DRB4** | **0.991** | | | | |

DRB4 needed one correction to get there: jMorp assigns every chromosome an allele and uses the
null `DRB4*03:01N` — 0.5430 of all chromosomes — for "this haplotype carries no functional
DRB4", which is what this pipeline records as `not_typed` and leaves out of its denominator.
That single entry was the entire discrepancy (`max|Δ|` = 0.5430 exactly); dropping nulls and
renormalising moves DRB4 from *r* = 0.625 to 0.991, and the carriage rates then agree
independently — ours 0.4018, jMorp non-null 0.4565.

**DRB3 is NOT COMPARABLE, and the reason is now known.** Earlier versions of this section said
three mechanisms had been tested and none accounted for the gap. That was wrong, and the
mechanism is the denominator.

jMorp carries no null allele for DRB3, yet its DRB3 counts still span all 122,848
chromosomes — every one in the panel. That is not possible: DRB3 exists only on DR52
haplotypes, and we call it on 47.7 % of control chromosomes. The panel is therefore assigning
a real DRB3 allele to chromosomes that do not carry the gene, and they land on the commonest
one: `DRB3*01:01P` is 0.674 in the reference against our 0.189. **Remove that one allele and
renormalise both sides and *r* moves from 0.635 to 0.944** — most of the gap is that single
absorbing allele.

So our denominator is "chromosomes carrying DRB3" and the reference's is "all chromosomes".
Those are different quantities and the correlation between them measures the mismatch, not the
typing. `allele_freq_check.py` now decides this **mechanically** — a locus is marked
`comparable = 0` when the reference covers essentially every chromosome while we call the gene
on materially fewer — and DRB3 is excluded from the headline and drawn hollow. DRB4 is not
caught by the same test and should not be: the 45.6 % of reference chromosomes it puts on
a real DRB4 allele is an honest fraction, not a saturated
one.

The remaining 0.944 is still below every comparable locus (lowest 0.983), so one mechanism
explains most of this and not all of it. What is settled is that the comparison at DRB3 was
never a measurement of typing quality.

### What replaced the 1000 Genomes panel, and why

The previous reference was **1000 Genomes JPT**: 105 individuals over A, B, C, DQB1 and DRB1.
At *n* = 105 the standard error near a frequency of 0.4 is ~0.034 — it bounds gross error and
nothing finer — and it carried no DPB1, no DQA1 and no DRB3/4/5, so **the loci with the least
reliable typing here were exactly the ones it could not check**. jMorp is ~580× larger and
covers eight of those previously-unchecked loci.

It is **kept and RUN**, not merely kept: `ALLELE_FREQ_CHECK` fans out over both panels on
every execution, each with its own locus list, writing `*.1kg.tsv` beside the primary
tables. Our calls reproduce it at *r* = 0.968–0.993 across the five loci it carries, against
*r* = 0.983–0.993 for the same five loci in jMorp. Two unrelated Japanese references
corroborating each other is worth more than either alone — and, as §3 shows, running both is
also the only way to demonstrate rather than assert that the unconfirmed *level* belongs to
the panel while the *contrast* does not.

**Neither is an accuracy measurement.** Both compare population frequencies; neither can say a
given sample was typed correctly, and a set of errors that happens to preserve the spectrum is
invisible to both.

### The measurement that would settle it — now the ONLY one available

Typing samples whose types are already published. Scoped on 2026-08-25 and deferred by
decision; **§3 has since retired the alternative.** The downsampling control recorded there
cannot work — cases and controls are at the same measured depth — so no within-data comparison
separates the technical gap from phenotype, and this is what is left:

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

## 3. Controls are typed measurably worse than cases, and it is not depth or the panel

**Found by the frequency check of §2 — this is the defect that check was built to detect,
and it is the largest open question here.** All numbers are for the `full_mainland`
cohort and are derived in `results/_run_info/facts.json`.

Per locus, the share of control chromosomes assigned to an allele the reference panel
never lists, beside the total shortfall across alleles the panel carries at frequency
≥ 0.05:

| locus | tail alleles | tail chromosomes | share | common-allele shortfall |
|---|---|---|---|---|
| DRB4 | 73 | 580 | **19.15 %** | 0.2083 |
| DRB3 | 67 | 424 | **16.69 %** | 0.1711 |
| B | 216 | 329 | 6.18 % | 0.0632 |
| A | 205 | 323 | 6.07 % | 0.0820 |
| G | 11 | 319 | 5.99 % | 0.0626 |
| C | 173 | 264 | 4.96 % | 0.0610 |
| E | 32 | 187 | 3.51 % | 0.0351 |
| DPB1 | 106 | 178 | 3.34 % | 0.0501 |
| DPA1 | 22 | 111 | 2.08 % | 0.0209 |
| DRB1 | 82 | 105 | 1.97 % | −0.0554 |
| F | 13 | 77 | 1.45 % | 0.0150 |
| DQB1 | 42 | 50 | 0.94 % | −0.0062 |
| DQA1 | 18 | 25 | 0.47 % | 0.0013 |

**The two columns match at the class I loci, and that is the finding there.** At A, B and
C the tail is not extra diversity sitting beside a correct common spectrum; it is taken
*out of* the common alleles, and the shortfall is close to the tail share.
**They do not match everywhere**, and the earlier version of this section over-generalised:
at DRB1 and DQB1 the shortfall is negative — our controls call those loci's common alleles
slightly *high* — so whatever moves the class I spectrum is not acting uniformly. DRB3 and
DRB4 are the largest shares in the table and are also the two loci §1 says are encoded
wrongly; read them as a symptom of §1 before reading them as a symptom of this.

### What it is not

- **Not the ambiguity resolution of §5.** Traced per allele on the earlier run: `B*44:302`
  appeared 67 times and exactly 1 of those came from an ambiguous call; `B*40:379`,
  `B*07:381`, `B*35:380`, `DRB1*14:253` and `A*24:50` came from ambiguity **zero** times.
  HLA-HD reports these as unambiguous calls.
- **Not the reference panel being too small to have seen the allele.** This is now tested
  rather than argued: the same calls are scored against two panels of very different size
  on every run. See "not the panel" below.
- **Not the samples not being Japanese.** The cohort is ancestry-restricted to
  `full_mainland`, so an allele missing from a Japanese panel can no longer be explained
  by the sample being from somewhere else. Before that restriction the same contrast was
  4.74 % against 3.57 %; after it, 4.63 % against 3.53 %. Removing the benign explanation
  left the finding where it was.
- **Not a parsing artefact.** The reference's own `*`-flagged and `/`-ambiguous values were
  a separate bug, fixed; it moved `n_ref_only` and left these numbers intact.

### What it probably is

HLA-HD choosing a nearly identical rare allele over the common one out of a 46,005-allele
dictionary — `B*44:302:01` instead of `B*44:03`, `B*40:379` instead of `B*40:02`. The
frequency table meant to break exactly these ties (`-f freq_data`) **is** passed, and it is
a global registry count rather than a Japanese one, which is the wrong prior for this
cohort.

### Why it is not fixed here

Every candidate fix changes genotypes, and the choice is a study decision, not a code one:

1. restrict the dictionary to alleles observed in East Asians, and re-type;
2. post-hoc collapse tail alleles onto their nearest common relative — cheap, reversible,
   and an invention nobody else uses;
3. accept it, and treat the 2-field allele test as the primary analysis only for alleles
   the reference panel confirms.

### IT IS DIFFERENTIAL BY PHENOTYPE, and that is the part that matters

Call rate is flat between the groups — which
is what made this look benign at first. **Typing accuracy is not.** The share of
chromosomes on a confirmed allele, per locus:

| locus | controls | cases | difference | Fisher *P* |
|---|---|---|---|---|
| DRB4 | 80.85 % | 90.91 % | +10.06 % | 1.84e-09 |
| DRB3 | 83.31 % | 90.78 % | +7.48 % | 0.000186 |
| A | 93.93 % | 98.86 % | +4.93 % | 4.08e-12 |
| B | 93.82 % | 97.95 % | +4.13 % | 5.17e-08 |
| C | 95.04 % | 99.09 % | +4.05 % | 3.68e-10 |
| DPB1 | 96.66 % | 99.66 % | +3.00 % | 6.73e-09 |
| E | 96.49 % | 98.18 % | +1.69 % | 0.00752 |
| DPA1 | 97.92 % | 99.32 % | +1.40 % | 0.0029 |
| DRB1 | 98.03 % | 99.20 % | +1.17 % | 0.0131 |
| G | 94.01 % | 94.87 % | +0.87 % | 0.352 |
| DQB1 | 99.06 % | 99.89 % | +0.83 % | 0.00741 |
| DQA1 | 99.53 % | 99.89 % | +0.36 % | 0.164 |
| F | 98.55 % | 97.04 % | -1.51 % | 0.00246 |
| **pooled** | **95.37 %** | **96.47 %** | **+1.10 %** | 1.89e-07 |

Twelve of thirteen loci point the same way. F is the one exception and is also the locus
where the reference is nearly monomorphic, so its denominator carries almost no
information.

### The cause is NOT depth

An earlier version of this section attributed the gap to sequencing depth, reading the
`15x` and `30x` **platform labels as if they were measured depths**. They are targets, and
the measured values do not agree with them:

| platform | label | measured depth (median) | unconfirmed share |
|---|---|---|---|
| NovaSeq | 30x | 30.93× | 1.32 % |
| DNBSeq-G400RS | 30x | 36.25× | 2.23 % |
| **DNBSeq-T7** (cases) | 30x | **18.43×** | **3.36 %** |
| **HiSeqX** (controls) | 15x | **18.72×** | **4.63 %** |
| DNBseq-G400RS | 15x | 14.17× | 10.66 % |

**Cases and controls are at the same measured depth.** Median 18.65× against 18.72×,
Mann–Whitney *P* = 0.194. And the gap survives matching on it: restricted to the 17–21×
window, where both platforms have a median near 18.5×, HiSeqX controls are at 4.48 %
unconfirmed against DNBSeq-T7 cases at 3.16 % — Fisher *P* = 3.58e-07 on 1,868 samples.

Depth *does* have a real effect, and it is visible where phenotype is held fixed. Within
cases only: DNBseq-G400RS at 14.17× gives 10.66 %, DNBSeq-T7 at 18.43× gives 3.36 %,
NovaSeq at 30.93× gives 1.32 %. **But that effect does not explain the case–control gap,
because there is no case–control depth difference for it to act on.**

### The cause is NOT the reference panel either, and this is now shown rather than argued

The obvious objection to everything above is that "unconfirmed" is defined by a finite
panel, so the whole measurement could be an artefact of what that panel happens to carry.
Both panels are therefore run on every execution, on identical calls:

| panel | individuals | loci | controls | cases | ratio |
|---|---|---|---|---|---|
| jMorp 61KJPN-HLA | 61,424 | 13 | 4.63 % | 3.53 % | 1.31 |
| 1000 Genomes JPT | 105 | 5 | 7.42 % | 5.40 % | 1.37 |

**The level moves by more than half again between the panels; the ratio barely moves.**
That is the whole argument for reading only the contrast, and it is now a measurement
rather than a caveat. Panel (d) of `figures/typing_confound.png` is this table.

Note which direction it runs. The smaller panel gives the *larger* unconfirmed share,
because 105 individuals have not sampled the rare Japanese repertoire. An allele absent
from 105 people is weak evidence; absent from 61,424 it is strong. The captions in this
component said the opposite for three months after the panel changed, which is what
`verify.sh` §7 and §8 now exist to prevent.

### Which output survives it

The table below was measured on the superseded 3,569-sample run and is **carried over
deliberately**: it compares the RESIDUES of one allele against another, which is a property
of the IMGT alignments and not of which samples were typed. Only the `pairs` column depends
on the cohort, and it changes only by which tail alleles happen to be observed.

Comparing each tail allele's residues against the commonest confirmed allele of its
1-field family:

| locus | pairs | median residues differing | share differing by ≤ 2 | positions |
|---|---|---|---|---|
| A | 396 | **2** | 51 % | 390 |
| C | 338 | **2** | 55 % | 417 |
| B | 395 | **5** | 41 % | 483 |
| DQB1 | 117 | 47 | 15 % | 274 |
| DRB1 | 144 | 85 | 22 % | 290 |

**Class I residues are largely insulated**: a mis-call to a rare allele of the same family
moves 2–5 of the locus's several hundred positions. **Class II residues are not** — DQB1 and
DRB1 tail calls move tens of positions — but their tails are also the smallest in the
per-locus table above, 0.94 % and 1.97 % of control chromosomes.

Order of robustness, most to least: class I residues → class II residues →
`allele_dosage.tsv`. Rare rows of `allele_dosage.tsv` should not be tested without checking
them against **`05.qc/allele_pgroup_map.tsv`** first — that is the crosswalk with the
`in_reference` flag, one row per allele the cohort carries.
`allele_frequency_check.tsv` is a per-locus comparison table and does not answer the
per-allele question, which is what [OUTPUTS.md](OUTPUTS.md) says about it.

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
