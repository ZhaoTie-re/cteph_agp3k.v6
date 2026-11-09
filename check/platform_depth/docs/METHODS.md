# Methods — `check/platform_depth`

Checks of every signal of the study against technical artifacts, each run in all three nested cohorts (narrow,
intermediate, full mainland; each a selection of the full-mainland samples): **read down-sampling** (does depth
create the signal?), a **platform-sensitivity analysis** (do the case platforms agree?), **SNP arrays** (does a
second technology in the same people agree?) and **control sources** (does one source carry the signal?);
descriptive follow-up counts (§9) stand behind the flags but are not a check. Every re-test uses the
main analysis' own model, variant or marker set and pre-specified threshold; only the genotypes (down-sampling)
or the sample subset (platform, array and control-source analyses) change.

Design decisions of the user: every case platform but DNBseq-G400RS 15x is down-sampled so that its median
effective depth equals the HiSeq X controls' median (2026-11-03; first, 2026-10-30, the deep platforms only),
with a baseline and three replicates for every platform (2026-11-03); platform strata share all controls; the
signals are re-tested, not the gene-test λGC; the genes are the Tier 1–2 genes (2026-10-30).

## 1. Signals

| Type | What is re-tested | Model and threshold (as the main analysis) |
|---|---|---|
| GWAS | the 168 signal variants of the 16 common-variant loci (16 leads, 3 alternative leads, 149 credible-set members; `integration/saige_signals/results/00.selection/variants.tsv`) | SAIGE 1.1.9 logistic mixed model (sex + 10 BBJ-mainland PCs, full GRM, LOCO), P = `p.value.NA`; P < 5 × 10⁻⁸ genome-wide, < 10⁻⁵ suggestive |
| Genes | STBD1, OSBP (Tier 1), RPH3A, JAK2 (Tier 2) in both variant sets (`analysis/assoc_gene/results/_comparison/robust_genes/tiers.tsv`); per cohort and set exactly the variants rvtest tested | rvtest CMC (+ cmcWald) and SKAT-O, sex + 10 PCs, mean imputation; Bonferroni over the cohort's gene scan, BH with the gene's new P substituted into the published scan |
| HLA | every tested allele and residue column (`hla.typing.v2/results/06.assoc_inputs`), headline DQB1:57:V; every residue position, headline DQB1 57 | alleles and residues: SAIGE step 2 on the published null model (mean imputation), GEC threshold per scan; positions: fixed-model omnibus LRT, Bonferroni over positions |

The effect is β = SAIGE's BETA with SE = 1/√var, so that β/SE reproduces the primary no-SPA P. SAIGE's own
SE column is not used: on rows where the saddle-point approximation was applied it equals
|BETA| / qnorm(1 − p.value), the one-sided quantile of the two-sided SPA P.

## 2. Samples, CRAMs, depth

- Sample table from the info sheet (`WGS_Platform` verbatim), cohort membership from each cohort's
  `pheno.tsv`, unrelated samples from the fixed-model `.fam`, rare-variant samples from tuning.rv's call set.
- CRAMs: `Cram_Path` followed to the real file; 72 full-mainland control links whose target directory was
  removed fall back to `/LARGE1/gr10478/pub/WGS/cram/Control_BBJ/` with indexes from hla.typing.
- **Effective depth** of every genome: `samtools coverage -q 20 -Q 10` (HaplotypeCaller's mapping- and
  base-quality filters), unmapped / secondary / QC-fail / duplicate reads excluded, over ten genome-wide 1 Mb
  probe windows (the July audit's) and over every calling interval. Read length from 1,000 primary reads.
- **Target** D* = the median probe depth of the full-mainland HiSeq X controls (16.99×). Every case platform but
  DNBseq-G400RS 15x is brought to it, so that each platform's median equals the controls' median:
  - each NovaSeq and DNBSeq-G400RS 30x case (the 89 deep cases) keeps the fraction f_i = D* / D_i of its read pairs;
  - every DNBSeq-T7 case (319) keeps one platform fraction, f = D* / (the full-mainland DNBSeq-T7 median, 17.60×)
    = 0.96528, which keeps the platform's own spread of depth (`scripts/t7_fractions.py`,
    `01.depth/fractions.all.tsv`);
  - DNBseq-G400RS 15x (median 13.77×, below D*; depth cannot be added) is not touched.

  The 408 down-sampled cases are full-mainland members; each cohort's re-test swaps the ones it holds (the 89 deep
  cases are in all three cohorts; DNBSeq-T7: narrow 299, intermediate 309, full 319).

## 3. Down-sampling arms

- **Arms:** baseline (f = 1, every read through the same new processing) and R = 3 replicates (seed r = 1, 2, 3)
  for every down-sampled case (user, 2026-11-03: one replicate count for every platform; the first run's ten
  replicates of the deep cases are kept in `results/_superseded/downsampling_R10/`). The depth effect is replicate
  vs baseline: identical processing, only the reads differ. Re-processing fidelity is baseline vs original.
- **Reads.** The calling intervals (§4) ± 1 kb are extracted once per genome (`samtools view --region-file`,
  a multi-region query that writes a read once). `samtools view --subsample f --subsample-seed r` keeps or
  drops a read pair by a hash of its name and the seed: mates stay together, and the same read names are kept
  in the variant regions and in the HLA reads. Each arm's sample is renamed `<id>__<arm>`.
- **Calling** is the cohort's discovery command: GATK 4.1.4.0 HaplotypeCaller `-ERC GVCF` with the cohort's
  settings (`--pcr-indel-model CONSERVATIVE`, MQ ≥ 20, BQ ≥ 10: the 4.1.x defaults its VCF header records).
  Alleles are not forced (`--alleles` turns reference blocks into site records and changes the GQ and DP of
  hom-ref genotypes, the main source of no-calls at control depth).
- **Joint genotyping**, per chromosome, all arms together with the surviving original GVCFs (one allele set):
  GenomicsDBImport (one partition per calling interval) and GenotypeGVCFs 4.1.4.0 twice — variant sites
  over the intervals with `-stand-call-conf 0`, and `--include-non-variant-sites` at the fixed positions
  (4.1.4.0 has no `--force-output-intervals`). `expand_fixed.py` writes one biallelic record per fixed site:
  the variant record of that allele; else 0/0 from another allele's record (AD = ref depth, DP and GQ kept,
  as `bcftools norm -m -any` gives a carrier of another allele); else 0/0 from the non-variant record
  (GQ ← RGQ, AD ← (DP, 0); without AD the AlleleFraction is NaN and the cohort's `AF == 'NaN'` rule would
  no-call a hom-ref genotype).
- **Genotype QC** is the cohort's: `bcftools norm -m -any --check-ref s`, drop `ALT="*"`, IDs
  `CHROM:POS:REF:ALT`; VariantAnnotator `-A AlleleFraction` and VariantFiltration from the cohort's GATK 4.6.2.0
  image; `+setGT` unphase; GQ < 20, DP < 8, het allele fraction outside 0.2–0.8, AF NaN → no-call. Site
  filters are not re-applied: the fixed sites passed them in the cohort.
- **Swap.** The deep cases' QC'd genotypes replace theirs in the original GWAS fileset (random-model
  variants) and gene fileset (tuning.rv call set), byte for byte (`bedio.py swap`); every other genome stays
  identical, and the variant set never changes.
- **HLA.** The archived HLA read pairs (`hla.typing/results/01.reads`, the reads that were typed: those on
  chr6:29.54–33.42 Mb, on the 16 chr6 ALT scaffolds and the 525 HLA decoy contigs of hs38DH, and the unmapped
  reads; `hla.typing/results/00.manifest/contig_list.txt`) are
  subsampled with the same fraction and seed (`samtools import` → `--subsample` → `samtools fastq`) and typed
  by HLA-HD 1.7.1 from a checksummed snapshot of the install that typed the cohort, with the same dictionary,
  frequency data and options. hla.typing.v2 steps 01 → 02 → 05 → 06 run unchanged on a side input tree
  (`hla_v2_wrapper.py`), and the re-test uses the ORIGINAL tested columns and thresholds (`hla_matrices.py`).

## 4. Calling intervals and fixed sites

Each signal variant (the 168 GWAS signal variants and the 36 tested variants of the four genes) ± 10 kb, merged
(23 regions, 930 kb; listed with their signals in the deck's appendix), clipped to the original calling intervals
(`wgs_calling_regions.hg38.<chr>.interval_list`) and aligned to their 5-kb read-shard grid (HaplotypeCaller
splits each interval into 5,000 bp shards from its start; on the grid the shards are the original run's).
Fixed sites are every variant the analyses read inside these intervals: random-model variants of any cohort
and the genes' tested variants.

## 5. Concordance (a reported result)

Genotypes are compared at every fixed site as ALT counts 0 / 1 / 2 or missing: baseline vs original,
original GVCFs re-genotyped by this procedure vs original (the joint-genotyping step), baseline vs those
(the HaplotypeCaller re-run), each replicate vs baseline (depth), each replicate vs original (total). Per
comparison: the 4 × 4 confusion counts, concordance among genotypes called in both, non-reference concordance,
het and hom-alt retention, loss to missing, false-het rate and call rates, by sample, platform, variant class
and MAF class; replicates as median and range. The share of genotypes each QC rule no-calls is compared
between the arms, the fidelity controls and the cohort's own pre-QC genotypes by group, in each cohort (median
over each group's genomes; `ds.qc_reasons.tsv`). HLA: per gene the
call state and the two-field allele copies shared, and the dosage of the tested columns, arm against archive
and baseline; the allele read support against the controls'.

## 6. Platform sensitivity

Subsets per cohort, every one with **all** the cohort's controls: the pooled set (every case; its files
equal assoc_hla.v2's 00.prep), four strata (one case platform each) and four leave-one-platform-out sets.
Each gets its own SAIGE null model (assoc_saige's settings and GRM markers). The pooled re-fit is checked
against the published null (gate G5).

- **(a) Heterogeneity.** With shared controls the stratum estimates are correlated:
  Σ_kl = ρ_kl SE_k SE_l, ρ_kl = √(n_k n_l / ((n_k + n₀)(n_l + n₀))) (Lin & Sullivan 2009; n_k cases of
  stratum k, n₀ controls), checked against the correlation of z over ~20,000 random LD-pruned SNPs and indels.
  GLS estimate β_FE = 1ᵀΣ⁻¹β / 1ᵀΣ⁻¹1 and Q = (β − β_FE)ᵀ Σ⁻¹ (β − β_FE) ~ χ²(K − 1) (Cochran 1954);
  I² = max(0, (Q − df)/Q) with the test-based CI of Higgins & Thompson (2002). Q's P is also read empirically
  against the null SNPs of the same MAF bin and type, because three of the four strata are small. Each stratum
  carries the power to detect a stratum with no effect.
- **(b) Meta-analysis vs pooled.** β_FE weights the strata by their case counts and so tracks the pooled
  estimate by construction (reported as a check). The random-effects estimate (τ² by REML with
  V = Σ + τ² I, Hartung–Knapp CI on t(K − 1)) is compared with the pooled result (Evangelou & Ioannidis 2013). Like Q,
  it is read at MAF ≥ 5 % only: in small strata SAIGE's one-step β of a rare variant is not reliable.
- **(c) Leave one platform out.** β, SE and P without each case platform, with the z expected if the effect
  were unchanged, β_pooled / SE_(−k): fewer cases lower the P without the effect changing. The difference
  β_k − β_(−k) is tested with the shared-control covariance; for NovaSeq it compares the only cases on the
  controls' chemistry (Illumina) with the rest (MGI).
- **Case-only check.** Within unrelated cases, multinomial logistic regression platform ~ marker + sex +
  10 PCs, LRT with K − 1 df: with shared controls, equal ORs across platforms means an equal marker
  distribution across the case platforms (Begg & Zhang 1994). Genes (0–2 carriers per stratum):
  carriers × platform, permutation χ²; per-stratum carrier ORs (Firth) are descriptive.
- Positions: the omnibus per stratum and leave-one-platform-out set; their heterogeneity is read through
  their residues.

An artifact shared by every case platform cannot be seen by these comparisons; the SNP arrays cover the
DNBSeq-T7 cases (§7), and the down-sampling covers the deep cases.

## 7. SNP arrays (`scripts/array_cohorts.py`, `scripts/array_evidence.py`, `11.array`)

The same people genotyped by a second technology. The unrelated genomes with SNP-array genotypes are DNBSeq-T7
cases (Illumina ASA24v1) and HiSeq X controls (JSA24v1); the NovaSeq and DNBSeq-G400RS cases have none:

| Cohort | People with arrays | Cases (DNBSeq-T7) | Controls (HiSeq X) |
|---|---|---|---|
| Narrow mainland | 1,398 | 279 | 1,119 |
| Intermediate mainland | 1,460 | 289 | 1,171 |
| Full mainland | 1,683 | 299 | 1,384 |

The arrays are confounded with case status too, but by another technology, so agreement argues against a WGS
artifact without excluding an error shared by both. `check/array_crosscheck` made the comparison in full
mainland; `scripts/array_cohorts.py` re-runs its steps on its inputs (read only) in each cohort, keeping that
cohort's unrelated members (fixed-model set) with its phenotype and covariates; nothing else changes (the shared
sites, the model, the procedures). Gate G10: in full mainland the re-run reproduces array_crosscheck's tables
(`11.array/array.gate.json`). `scripts/array_evidence.py` then files the result per signal and cohort, comparing
WGS with the array in the same people, at the same marker and allele.

- **Genome-wide.** Genotype concordance of non-missing calls at the shared sites, cases and controls apart; the
  same logistic model (sex + 10 BBJ-mainland PCs, Firth fallback; plink2) on both technologies at the QC'd
  shared sites (MAF ≥ 0.01, HWE P ≥ 10⁻⁶ in controls on both): λ of each and the correlation of z.
- **One marker per signal** (user, 2026-11-05: a simple rule; cases and controls are read from the same kind of
  data). The lead typed on the arrays; else the lead in the imputed arrays (Minimac4 v4.1.4 against the AGP3K
  whole-genome panel, chromosomes 1–22) with imputation R² ≥ 0.8 (Minimac's own R², one value per marker; the
  imputation ran once, so it is the same in every cohort); else a typed LD proxy with r² ≥ 0.8 in the cohort's
  controls; else not checkable. DQB1:57:V: DQB1 imputed with HIBAG (`scripts/array_hibag.py`): HIBAG 1.42.0 with its
  published Asian two-field model (`assets/hibag/`, md5 in `MD5SUMS`; 612 training people of East and South Asian
  ancestry including HapMap CHB and JPT, 500 classifiers; `DQB1_model.json`), on the QC'd array SNPs of the xMHC
  (chr6:28–34 Mb) called in ≥ 98% of the cases and of the controls; residue 57 from hla.typing.v2's IMGT reference.
  Each marker carries one Minimac R² and its genotype concordance with WGS in cases and in controls apart (user,
  2026-11-05): the share of people called on both whose array genotype equals the WGS genotype (typed: the array
  call; imputed: the best guess, dosage < 0.5 → 0, < 1.5 → 1, else 2; DQB1:57:V: HIBAG's best-guess 57 V count;
  `06.imputed/concordance.tsv` also gives it among people non-reference on either side). The controls are in the
  imputation panel, so their concordance is near 1 by construction. Odds ratios are given per the lead's ALT allele
  (array_crosscheck counts its A1, often the minor allele; the deck reads every OR per ALT).
  The association is re-run on that marker in the same people (sex + 10 PCs) and compared with WGS by z (signed by
  the OR, sized by the two-sided P): contradicted = opposite direction or |z_WGS − z_array| > 3; otherwise
  consistent (user, 2026-11-05: an array P < 0.05 is not required, since the arrays hold only the DNBSeq-T7 cases and
  power alone can leave a real signal above 0.05; the array P is kept in the tables). The imputation panel holds the controls' own genomes, so their imputed genotypes largely repeat WGS;
  the check rests on the cases.
- **Gene sets** (user, 2026-11-05: a description, no verdict, no expected number). The arrays cannot reproduce a
  rare-variant gene test: the people are the same, the controls are in the imputation panel, and imputation misses
  many rare carriers. The drivers are the six of the gene-based analysis (the leading variant of each tiered M/H set
  and any variant whose removal lowers −log₁₀ P by ≥ 1 in all six tests; `array_evidence.py` checks them against
  `integration/gene_signals/results/gene_signals.lovo.tsv`); array_crosscheck's own list also held STBD1 p.Ser188Arg,
  which is not one (its removal raises −log₁₀ P), and stays in its tables (gate G10) but not in the readout.
  What they can show is how many WGS case carriers of each driver the imputed arrays confirm (dosage ≥ 0.5): the
  cases' imputed genotypes come from their own arrays and the panel haplotypes only, without their sequencing reads,
  so a confirmed carrier is evidence, independent of the case's own reads, that the WGS call is real; a missed one is
  weak evidence against,
  since rare-variant imputation misses carriers in general (41–83% of the case carriers of the rare untyped variants
  within 500 kb of the drivers are recovered, by their number of control carriers, `06.imputed/recovery_by_group.tsv`;
  a context number, not a test). A case imputed as a carrier without a WGS call of it is counted too. Controls are
  in the panel and carry no information (99.6% of their carriers of the same rare variants are recovered). A driver
  absent from the imputed arrays (JAK2 p.Val617Phe: typed on both array products as a multi-allelic site and removed
  at array QC; most likely acquired, absent from the panel) is not checkable.
  The carriers are counted once, in full mainland, which contains the narrow and intermediate cohorts; the matrix
  shows one 'reported' cell for the three cohorts. `array.evidence.tsv` keeps an expected count per driver (WGS case
  carriers × the recovery at its number of control carriers) as a diagnostic; it is not part of the readout.
- **Per signal.** Contradicted if any row is; array_cohorts.py still re-runs array_crosscheck's tables (including
  its ridge prediction of 57 V) for gate G10 only.

It is evidence beside the platform and depth checks, limited to the DNBSeq-T7 cases and the HiSeq X controls.

## 8. Control sources

The controls come from five projects (the info sheet's Project column): Nagahama (a population cohort), BBJ,
HTLV-1, ACC and lung cancer; random-model controls per cohort:

| Cohort | Nagahama | BBJ | HTLV-1 | ACC | Lung cancer | All |
|---|---|---|---|---|---|---|
| Narrow mainland | 1,075 | 195 | 168 | 198 | 140 | 1,776 |
| Intermediate mainland | 1,108 | 342 | 241 | 227 | 161 | 2,079 |
| Full mainland | 1,298 | 634 | 310 | 242 | 178 | 2,662 |
 Each cohort gets ten more subsets,
rows of the pooled platform subset (`scripts/control_subsets.py`): every case with one source's controls
(`ctrlonly_<source>`) and with all sources but one (`ctrlwithout_<source>`). They run through the platform
machinery, each signal with its main model: an own SAIGE null model and step 2 (GWAS signal variants, HLA alleles
and residues), rvtest CMC and SKAT-O (Tier 1–2 genes, both variant sets), the fixed omnibus (positions).
`scripts/control_report.py` reads the published results:

- **Across the five sources.** GLS, Q (its analytic χ² P is read, at MAF ≥ 5 %; over the null SNPs it is
  lenient, below 0.05 more often than 5 %, so an empirical P as in §6 would be larger; review 2026-11-06. The deck
  quantifies it from the subsets' SAIGE rows with this section's formulas (`deck/collect_stats.py`, `sources_null_q`,
  `sources_q_emp`): at MAF ≥ 5 % the analytic P is below 0.05 in 8.1 / 9.1 / 10.9 % of the null SNPs (N / I / F), and
  read against the null SNPs of its MAF bin and type the smallest source Q P, ABI3BP's 0.022 / 0.020 / 0.049, is
  0.042 / 0.047 / 0.107; the rule keeps the analytic P, which can only add flags), I², REML with Hartung–Knapp as in §6;
  the source estimates share every case,
  ρ_kl = √(m_k m_l / ((m_k + n)(m_l + n))) with m a source's controls and n the cases (Lin & Sullivan 2009 with
  the roles of cases and controls swapped).
- **Leaving each source out.** β, SE, P, z observed and expected (β_pooled / SE without the source), and the
  difference between that source alone and the rest.
- **Calibration.** λ and the share of |z| > 3 over the null SNPs in every subset.
- **Genes.** P, CMC OR and carriers per subset; per source the carriers are descriptive (0–3 per source).

For a rare variant in a small source subset SAIGE's one-step β is not reliable (§6): leaving one source out is
read first.

## 9. Follow-up counts (`scripts/signal_followup.py`, `13.followup`)

Counts and frequencies in each of the three cohorts, read behind the flags; they are descriptive, not a check of the
evidence matrix (2026-11-03), and carry no new association test (no per-variant test of a rare variant).

- **GWAS leads by group.** The ALT frequency in each case platform, split by the info sheet's `Registration_Study`
  (new / old), and in each control source, beside ToMMo 61KJPN and 1KG JPT. A platform-specific frequency whose
  genotypes the arrays confirm points to the case group, not to the calls.
- **Gene carriers and carriers set missing.** For every tested variant of the Tier 1–2 genes and for each gene:
  carriers per case platform and per control source. The cohort's genotype filters (GQ < 20, DP < 8, heterozygote
  allele fraction outside 0.2–0.8, allele fraction missing → no-call; `wgs.auto.par/select.auto.par.v6.nf:252-260`)
  are re-applied to the pre-QC records (`wgs.auto.par/results/03_annotate_af_norm_gt`) and must reproduce the
  published genotypes in every cohort (`followup.gene_filters.tsv`; checked by `verify.sh`). A genotype set missing although its record is
  non-reference or carries ≥ 2 ALT reads is a possible carrier hidden as missing: their number bounds how much of a
  case–control difference missingness can make (CMC counts a missing genotype as a non-carrier). Reference
  frequencies: ToMMo 61KJPN, gnomAD v4.1 EAS, 1KG EAS and JPT (`integration/gene_signals/results/03.frequencies`),
  and the carriers ToMMo predicts (Σ 2·AF).
- **The GWAS leads before and after the filters.** At every lead, per group: carriers in the pre-QC records and in the
  published genotypes, genotypes set missing by the filters, and those set missing although the record carries the
  ALT (`followup.lead_preqc.tsv`); and every pre-QC record within 50 bp of the lead, with its carriers per group
  (`followup.lead_window.tsv`): another representation of the lead's allele, or a linked variant, that a platform
  calls instead.
- **DQB1:57:V across the replicates.** Every down-sampled case whose 57 V dosage differs between any two arms, with its DQB1
  call, read support and HLA-HD confidence: a systematic depth effect moves many samples one way, a typing error
  moves one sample in some replicates. The deep cases belong to all three cohorts, and their typings are the same
  in each.

## 10. Evidence matrix (`scripts/evidence_matrix.py`, `14.evidence`)

Every signal (the 16 GWAS leads, DQB1:57:V, position DQB1 57, the Tier 1–2 genes in their M/H sets) × every
check × cohort (user, 2026-11-05: the M/H sets only, the main analysis; the L/M/H sets are re-tested like the M/H ones
and their tables are kept as data, but no reading, flag or count is made from them): the metric and a reading by a fixed rule (`evidence.rules.tsv`; proposed, not yet confirmed by the
user). The checks are kept once, in the registry `scripts/checks.py`: per check the artifact it targets, the samples,
the readout, the rule and its thresholds, the tables and the figure ([PIPELINE.md](PIPELINE.md)); the matrix, the
figure and `verify.sh` read their rules, thresholds and vocabulary from it. Readings: consistent; flag; not testable (no left-out group with |z expected| ≥ 2.8); not applicable (Q at
MAF < 5 %; per-stratum estimates of genes, whose strata hold 0–7 carriers). The array check reads consistent or not
checkable, and flags a contradiction (opposite direction or |z_WGS − z_array| > 3). Every check is read in each of the three
cohorts and flags a signal when its rule holds in any one of them (was two of three before 2026-11-04). **The matrix flags and never removes**: which flags keep a locus out of the paper is the
user's decision (2026-10-31: no signal is excluded for now).

The leave-one-out rule takes the pooled β as the truth: z expected = β_pooled / SE without the group, and a group
whose |z expected| ≥ 2.8 leaves |z observed| < 1.96 or the opposite sign is a flag (the rule as pre-specified).
Its properties: the estimate without a group is nested in the pooled one, so if the effect is the same in every
group, z without the group, given the pooled estimate, is normal with mean z expected and SD
√(1 − (z expected / z pooled)²). The winner's curse inflates both estimates alike and does not bias the comparison
(Zhong & Prentice 2008 and Xiao & Boehnke 2009 concern the pooled estimate itself). On the 20,000 null SNPs the
model holds for every platform left out; for the control sources it is approximate (each subset has its own null
model, and the SE without a small source can fall below the pooled SE), so the per-source chances are nominal
(review 2026-11-06; quantified in `deck/collect_stats.py`, `loo_null_sd`: the SD of (z − z expected) / SD over the
null SNPs is 0.93–1.01 for the platforms left out and 1.10–1.92 for the sources, and the SE without a source is not
above the pooled SE in up to 52 % of the null SNPs, Lung cancer in full mainland). The SD is about 0.8 when a large
group (DNBSeq-T7) is left out and about 0.3 for a small one, so the rule is sensitive where a large group is left out
but such a group often fails the testability bar (its z expected is low), while a small group always passes the bar
and moves z little: an excess confined to a small group is described by the difference test (`p_diff`); for a rare lead in a small
group `p_diff` rests on the one-step stratum β, which is unreliable (§6), so the allele counts are read instead (the
deck adds Fisher's exact test of one platform's ALT alleles against the other cases', from `followup.gwas_af.tsv`, for
every rare lead and case platform: `deck/collect_stats.py`, `rare_platform_exact`; descriptive, e.g. INPP5D in the
DNBSeq-G400RS 30x cases P 0.019 in full mainland, where `p_diff` is 1.6 × 10⁻⁶). The number
of flags expected by chance under this model is reported beside the observed flags
(`evidence.loo_calibration.tsv`). The difference test and the follow-up counts (§9) are read with every flag.

Why the gene sets are read with the CMC burden β, and why JAK2 is not testable (user, 2026-11-06: a rationale for the
rule as it is). (1) z expected needs an effect and its SE; a burden test gives one (the CMC Wald β), SKAT gives none
(it tests the variance of the variant effects, with no direction) and SKAT-O only through its burden weight ρ. (2) The
test was fixed for every set before the readings; choosing it per set after the results would be post hoc. (3) The
P without a group is read, but only where it can be interpreted: |z| < 1.96 is P > 0.05 without the group, and the rule
reads it only where |z expected| ≥ 2.8, so that an unchanged effect would keep P < 0.05 at least 80 % of the time;
elsewhere the rule does not read a lost P. Given the pooled estimate the spread is smaller than 1 (0.20–0.86), so
the bar is conservative: three readings below it would keep P < 0.05 with 81–91 % probability (HS3ST3A1 without
DNBSeq-T7 in I and F; RPH3A M/H without Nagahama in N); the deck discloses them with their chance, and the bar, fixed
before the readings, is kept (review 2026-11-06). A SKAT-O P without a group has no expected value under an unchanged effect, so it
cannot be read this way. (4) Hence JAK2 M/H, called by SKAT-O with ρ = 0 (pure SKAT) and without a burden effect (CMC P
0.17 / 0.34 / 0.22, N / I / F), is 'not testable', not 'consistent'. Its SKAT-O P with each platform and each control
source left out (`09.platform_assoc/<cohort>/gene.*.moderate_high.tsv`) and its case carriers by platform are
described beside the rule (deck): strongest without DNBSeq-T7, which holds 4 of its 10 case carriers (1.3 % of its
cases against 5.0 % elsewhere). (5) The depth check reads each set with its calling test (SKAT-O for JAK2 and RPH3A):
it compares the same people before and after, so a P needs no expectation.

## 11. Gates

No comparison is read before the machinery reproduces the published results (results/_gates):
G1 SAIGE step 2 on the original genotypes reproduces `saige.assoc.tsv`; G2 rvtest reproduces the published
gene P, OR and carriers; G3 the hla.typing.v2 wrapper reproduces the six 06 matrices (md5); G4 the HLA scripts
reproduce `assoc_all.tsv` and the omnibus; G5 the pooled null re-fit reproduces the variance ratio and z
(1 %); G6 re-processing fidelity; G7 the swap leaves every genome but the deep cases' byte-identical
(`bedio.py swap` aborts otherwise; `04.genotypes/<label>/<arm>/<cohort>/swap.*.json`, checked by `verify.sh`); G8
achieved depth; G9 HLA re-typing at f = 1; G10 the array comparison re-run in full mainland reproduces
check/array_crosscheck's tables. Rules and records of every gate: [PIPELINE.md](PIPELINE.md#gates).

## References

Begg CB, Zhang ZF (1994) Cancer Epidemiol Biomarkers Prev 3:173–175 · Cochran WG (1954) Biometrics 10:101–129 ·
Evangelou E, Ioannidis JPA (2013) Nat Rev Genet 14:379–389 · Garner C (2011) Genet Epidemiol 35:261–268 ·
Hartung J, Knapp G (2001) Stat Med 20:3875–3889 · Higgins JPT, Thompson SG (2002) Stat Med 21:1539–1558 ·
Lin DY, Sullivan PF (2009) Am J Hum Genet 85:862–872 · Xiao R, Boehnke M (2009) Genet Epidemiol 33:453–462 ·
Zhong H, Prentice RL (2008) Biostatistics 9:621–634.
