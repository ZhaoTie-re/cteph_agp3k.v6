# hla.typing.v2 — outputs

All tables are tab-separated with a header row. "Reported genes" are the 19 functional genes
(*A, B, C, E, F, G, DRA, DRB1, DRB3, DRB4, DRB5, DQA1, DQB1, DPA1, DPB1, DMA, DMB, DOA, DOB*). The 14
pseudogenes HLA-HD also types appear in `calls_long.tsv` only.

## 00.inputs — the archived typing run

| file | content |
|---|---|
| `typing/<sample>.tar.gz` | 75 files: `result/` (final result, 33 × `est.txt`, 33 × `read.txt`), `hlahd.run.log`, `typing_validation.tsv`, `task/` (HLA-HD and extraction scripts and logs) |
| `archive_manifest.tsv` | sample_id, md5, n_members, typing_source, extract_source (paths into `hla.typing/work/`) |
| `sample_inputs.tsv` | sample_id, hlahd_version, read_pairs, regions, singletons_discarded, reads_processed, validation, group, sex, platform, observed_depth, archive_md5 |
| `provenance.json` | version, the two commands, read-pair summary, md5 of every pinned file |
| `hlahd_freq_data/` | HLA-HD's allele-frequency table, as used |

## 01.calls — genotypes

**`calls_long.tsv`.** One row per sample × gene (3,101 × 33).

| column | meaning |
|---|---|
| reported | 1 for the 19 reported genes |
| state | heterozygous, homozygous, ambiguous, unresolved, failed; for *DRB3/4/5*: present, absent, missing, ambiguous |
| copies | 2 for obligate genes; 0, 1 or 2 for *DRB3/4/5* |
| allele_1, allele_2 | the 2-field alleles (blank where none); for one-copy *DRB3/4/5*, allele_2 is blank |
| basis_1, basis_2 | agreed, prior, pgroup or none (METHODS §3.3) |
| p_group_1, p_group_2 | P-group representative |
| agreed_1, agreed_depth_1, agreed_2, agreed_depth_2 | the deepest name all candidates share, and its number of fields |
| raw, raw_result | HLA-HD's call type (one, same, two, not_typed, missing) and its result line |
| homozygous_model | 1 if HLA-HD's best pair was a single allele |
| n_best, n_other, n_shown | tied best pairs; further ambiguous pairs; entries on the result line |
| n_cand_1, n_cand_2 | candidate allele names per side |
| prior_ratio_1, prior_ratio_2 | for basis = prior, the frequency table's count ratio (chosen / best other) |
| reads_1, reads_2, unique_1, unique_2, shared | reads assigned to each allele, unique to each, and shared |
| drb_expected, drb_observed, drb_flag, raw_allele_1, raw_allele_2 | *DRB3/4/5* reconciliation (METHODS §3.5) |

**The other tables.**

| file | content |
|---|---|
| `genotypes.tsv` | wide: one row per sample, `<gene>_1`, `<gene>_2`, `<gene>_state` for the reported genes |
| `tie_candidates.tsv` | every call with tied best pairs: the pairs, the reported pair, the basis |
| `drb_reconciliation.tsv` | per sample and paralogue: *DRB1* alleles, expected and observed copies, HLA-HD's result, final state and alleles, flag |

## 02.evidence — read evidence and confidence

| file | content |
|---|---|
| `call_evidence.tsv` | reported genes: state, alleles, basis, prior ratios, reads, unique reads, min_allele_reads, unique_total, minor_share, balance_assessed, low_support, imbalanced, drb_flag, weaknesses, confidence |
| `gene_summary.tsv` | per gene: states, 2-field sources, low_support, imbalanced, confidence classes, share_high |
| `sample_qc.tsv` | per sample: group, sex, platform, observed_depth, read_pairs, failed_obligate, no_2field, low_conf, high_conf, flags, outlier, outlier_reason |
| `drb_consistency.tsv` | gene × expected × observed copies, before reconciliation |

## 03.panels — external panels (controls only)

**`locus_concordance.tsv`.** One row per panel and locus.

| column | meaning |
|---|---|
| n_chr_ours, n_chr_panel, n_alleles, n_shared, r | chromosome counts, allele counts, and frequency correlation |
| D, D_null_mean, D_null_p95, excess_D, p_D | discordance against the sampling expectation |
| D_shared, D_shared_high | discordance over the shared alleles: all calls and high-confidence calls |
| unconfirmed, unconfirmed_assessable, unconfirmed_after_vintage, unconfirmed_high, unconfirmed_assessable_high | absent-allele shares (METHODS §5) |
| panel_vintage | newest IPD-IMGT release among the panel's alleles |
| frequency_comparable, note | 0 and the reason for ToMMo *DRB3* and *DRB4* (`config.TOMMO_FREQ_NOT_COMPARABLE`) |

**The other tables.**

| file | content |
|---|---|
| `allele_concordance.tsv` | panel, gene, P group: counts and frequencies on both sides, difference, Fisher P, BH q |
| `homozygosity.tsv` | obligate genes: observed and expected homozygosity, ratio and CI, dropout and CI, permutation P, ToMMo-based expectation |
| `panel_ceiling.tsv` | D, excess D and r: ToMMo vs 1KG JPT, and ours vs 1KG JPT |
| `shift_direction.tsv` | alleles ≥ 1 point apart from ToMMo at the 1KG loci: the three frequencies, and which side 1KG JPT is on |
| `unconfirmed_alleles.tsv` | every called allele absent from ToMMo: P group, count, count among high-confidence calls, first IPD-IMGT release, assessable |
| `confidence_vs_panel.tsv` | absent-allele share by call confidence, pooled and per gene (assessable alleles) |
| `excluded.tsv` | controls left out of a locus comparison, with the reason |

## 04.confidence

`locus_tiers.tsv`. Per reported gene:
- the criteria values: unconfirmed (assessable), unconfirmed_high, dropout and CI, share_high;
- descriptive columns: excess_D, D_shared, excess_D_1kg, unconfirmed_any, unconfirmed_after_vintage;
- pass flags (blank = not applicable), criteria_failed, tier.

## 05.exports — for `analysis/assoc_hla` (not assessed here)

| file | content |
|---|---|
| `allele_dosage.tsv` | sample × 2-field allele, 0/1/2; blank = no genotype at that gene |
| `residue_dosage.tsv` | sample × `gene:IMGT position:residue`, 0/1/2; polymorphic positions only; blank where the sample has no genotype or the alignment leaves one of its alleles undetermined at that position |
| `allele_pgroup_map.tsv` | gene, allele, p_group (ToMMo-style name), in_reference, freq_reference (ToMMo, nulls dropped) |
| `call_confidence.tsv` | sample × gene confidence class |
| `residue_sites.tsv`, `allele_match_map.tsv`, `allele_unmatched.tsv`, `residue_reference/`, `residue_reference_summary.tsv` | residue matching record |

## figures

See `results/figures/README.md`.

## 06.assoc_inputs — what the association reads (no test is run here)

| file | content |
|---|---|
| `<cohort>/allele_dosage.tsv` | the cohort's samples (phenotype-file order) × allele columns that pass every rule; 0/1/2, blank = no genotype |
| `<cohort>/residue_dosage.tsv` | the same for residue columns |
| `markers.tsv` | one row per cohort and exported column: class, gene, position, residue, samples, called, call_rate, call_rate_cases, call_rate_controls, copies, mac, maf, status (the call-rate rule uses the two group columns) (`not_varying` / `excluded_gene` / `call_rate` / `maf` / `tested`) |
| `qc_record.tsv` | columns left after each step (`exported`, `genes`, `call_rate`, `maf`), per cohort, gene and class |
| `summary.json` | thresholds; per cohort the samples, cases, controls, columns tested and the md5 of each matrix |

Cohorts: `narrow_mainland`, `intermediate_mainland`, `full_mainland`. Columns tested (alleles + residues):
143 + 672, 144 + 681, 142 + 679.
