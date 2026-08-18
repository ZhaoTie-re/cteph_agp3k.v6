# Study notes — this dataset's run of `assoc_plink2`

**`assoc_plink2` is a reusable component.** [`METHODS.md`](METHODS.md) describes the *method* and holds
nothing specific to any dataset. This file is the other half: every finding, measured value, named
locus and configured resource that belongs to **this** project's run. Nothing here is required to run
the component elsewhere, and nothing here should be copied into a different study.

---

## 1. Configuration used

| parameter | value |
|---|---|
| `Cohorts` | `narrow_mainland`, `intermediate_mainland`, `full_mainland` (nested, narrowest first) |
| `Models` | additive, dominant, recessive; `PeakModel = additive` |
| `PcLabel` / `NPcs` | `bbj_mainland` / 10 — projections onto a mainland-restricted BBJ PCA space |
| `CovarLabel` | `SEX + 10 bbj_mainland PCs` (11 terms) |
| `FirthMode` | `no-firth` |
| thresholds | `PGenomeWide = 5e-8`, `PSuggestive = 1e-5`, `PeakFlank = 250 kb` |
| phenotype | CTEPH case/control |

The PCs describe ancestry structure only. They are **not** sequencing-batch or depth covariates.

### External resources

| purpose | path | version |
|---|---|---|
| PLINK 1.9 (`--r`/`--r2`) | `/home/b/b37974/plink` | v1.9.0-b.7.7 |
| PLINK 2 (GLM, geno-counts, hardy, window extraction) | `/home/b/b37974/plink2` | v2.0.0-a.5.17LM |
| tabix | `/home/b/b37974/htslib-1.9/tabix` | htslib 1.9 |
| rsID VCF (`rsid_vcf`) | `ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz` | 60KJPN |
| population co-occurrence r² (`PopR2Template`) | `ToMMo_60KJPN/co-occurrence/tommo-54kjpn-*-chr{N}-plink-r2.tsv.gz` | 54KJPN |
| reference panel (`RefPanelBfile`) | `cteph_agp3k/.../06.regional_plot_rev1/eas_all` | 1000G-EAS, 67.4 M × 504 |
| recombination track | `…/06.regional_plot_rev1/info/recomb1000GAvg.bw` | 1000G average |
| gene models | `EnsDb.Hsapiens.v86` via conda `r_work` | GRCh38, Ensembl 86 |
| snpEff index | `JHRPv6/.../snpEff.v6.index` | platform pre-computed |

`LdPanelLabels` name these on the regional figure: `ToMMo 54KJPN`, `1000 Genomes EAS (n = 504)`.

**No dbSNP build exists on this system**, and the reference panel and every cohort/source VCF carry `.`
in the ID column, so the ToMMo AF VCF is the only rsID source available. It is autosome-only, which is
not a limitation for an autosomal scan.

## 2. Scan size and calibration

| cohort | cases | controls | N_eff | λ_ADD | λ_DOM | λ_REC | variants analysed |
|---|---|---|---|---|---|---|---|
| `narrow_mainland` | 419 | 1,749 | 1,352 | 1.120 | 1.059 | 0.670 | 5,113,526 |
| `intermediate_mainland` | 429 | 2,051 | 1,419 | 1.121 | 1.060 | 0.665 | 5,113,902 |
| `full_mainland` | 439 | 2,632 | 1,505 | 1.132 | 1.065 | 0.654 | 5,115,079 |

### λ_ADD ≈ 1.12 is not noise, and it is not corrected

Relateds are already removed from the fixed-model set and 10 global PCs are already fitted, so what
remains is fine-scale structure that global PCs do not absorb — and case status is correlated with it.
The PopGMM keep lists quantify why:

| set | GMM components kept | cases kept | controls kept |
|---|---|---|---|
| `narrow_mainland` | 9 of 17 | **95.4 %** | **66.7 %** |
| `intermediate_mainland` | 12 of 17 | 97.7 % | 78.1 % |
| `full_mainland` | 17 of 17 | 100 % | 100 % |

Nine of seventeen mainland components retain 95.4 % of cases but only 66.7 % of controls — exactly the
condition a **GRM random-effect model** is designed to absorb and a fixed set of global PCs is not.

Consequences: the peak list is a **candidate** list; λ is reported and never applied; the follow-up is
the `random_model` genotypes (already built at `12_model_inputs/<cohort>/genotype/random_model/`)
analysed with **SAIGE** and **REGENIE**, as a separate component.

**LDSC is deliberately NOT added.** It is the standard way to split confounding from polygenicity, but
it needs samples in the tens of thousands; at N_eff ≈ 1,400 the intercept's standard error swamps the
estimate. Recorded so it is not added later by reflex.

### λ_REC ≈ 0.65 is deflation, not calibration

The recessive coding leaves most variants with too few alt-homozygotes to fit; plink2 drops 6.1–8.4 %
of the call set on `VIF_INFINITE` (310–440 k variants) and what survives is tested against a null that
is too narrow. Recessive peaks are read against this, not beside the additive results.

### Sequencing platform is confounded with phenotype

All cases are DNBSeq/NovaSeq, all controls HiSeqX 15×, so platform cannot enter as a covariate. Three
checks exonerate it as the driver of the inflation: no directional AF bias against the population
reference; case-vs-control AF-difference inflation 1.11 against 1.62 for a same-platform comparison;
and the risk-direction skew confined to low MAF.

## 3. Peaks called

| cohort | additive genome-wide | additive suggestive |
|---|---|---|
| `narrow_mainland` | 0 | 40 |
| `intermediate_mainland` | 1 | 37 |
| `full_mainland` | 2 | 51 |

Over M ≈ 5.11 M analysed variants a 1 × 10⁻⁵ threshold yields ≈ 51 variants by chance alone; the
observed counts run ~4× that, consistent with λ_ADD ≈ 1.12 rather than with 40–50 real loci. This is
why the suggestive tier receives no per-peak follow-up.

### Additive genome-wide loci

| locus | lead | cohort | *P* | OR (95 % CI) |
|---|---|---|---|---|
| *FTO* | `chr16:53887925:T:C` / rs16952623 | `intermediate_mainland` | 1.15 × 10⁻⁸ | 1.648 (1.388–1.957) |
| *FTO* | same | `full_mainland` | 2.88 × 10⁻⁸ | 1.591 (1.350–1.874) |
| *FTO* | same | `narrow_mainland` | 7.39 × 10⁻⁸ — **suggestive** | 1.631 (1.365–1.950) |
| *MIR3681HG* | `chr2:11891017` | `full_mainland` | 4.34 × 10⁻⁸ | — |

Credible sets: *FTO* 6 variants in both `intermediate_mainland` and `full_mainland`; *MIR3681HG* 21.

### Non-additive loci (reported, no fan-out)

- `full_mainland` REC: 7 genome-wide variants merging into **one** locus, `chr2:144936354`
  (rs11682379, *TEX41*, intronic, MODIFIER). Its lead departs from HWE in cases
  (`Case_HWE_P` = 5.2 × 10⁻⁴) while controls do not (0.78) — an excess of alt-homozygotes in one
  group, which is what makes a recessive test fire.
- `intermediate_mainland` DOM: rs4700094, *SREK1*, OR 0.52 (protective), *P* = 4.6 × 10⁻⁸.
- `intermediate_mainland` REC: rs12587119, *GALC*.

## 4. Defects this data exposed

These are the worked examples behind rules that [`METHODS.md`](METHODS.md) states abstractly.

**Degenerate fits with a clean ERRCODE** (METHODS §5). At a near-fixed variant the dominant coding
gives complete separation and plink2 returns `OR = inf`, `P = 0`, `ERRCODE = .`:

```
chr2:24542268:T:TAA   A1_FREQ 0.989   OR=inf   SE=4.19e6   Z=1.14e6   P=0   ERRCODE=.
chr10:54207183:A:G    A1_FREQ 0.980   OR=inf   SE=4.19e6   Z=1.61e7   P=0   ERRCODE=.
chr20:11681469:G:A    A1_FREQ 0.986   OR=inf   SE=2.97e6   Z=4.21e6   P=0   ERRCODE=.
```

Because `0 < 5e-8` each satisfied any significance test. The same three are entirely null additively
(*P* = 0.62, 0.85, 0.84). Incidence: 3 in `full_mainland/dominant`, 1 in
`intermediate_mainland/dominant`, 0 elsewhere. λ_GC was never affected — that computation already
required *P* > 0.

**One variant, two models, two effects** (METHODS §8b). `chr5:66174098` is a dominant genome-wide peak
(*P* = 4.6 × 10⁻⁸, OR 0.52) *and* an additive suggestive one (*P* = 1.7 × 10⁻⁷, OR 0.59). Annotating a
union keyed by variant ID handed one model's effect to the other, and really did print the additive
*FTO* peak with the dominant odds ratio (1.684 against 1.591).

**Blank ≠ not tested** (METHODS §12). The earlier comparison figure dropped rows where a cohort had not
called a peak, so *MIR3681HG* appeared with two rows instead of three.

**Truncated LD read as missing** (METHODS §10c). Verified in a 20 kb window at the *FTO* locus: 35,271
pairs, minimum *r²* = 0.200054, none below — the floor is a bound, not missing data. At the *FTO* peak
the co-occurrence panel is 1,112/1,112 informative, of which 1,053 are bounded below 0.2; the earlier
implementation rendered all of those grey and made the source look useless.

**Gene track drowned by clone accessions** (METHODS §9). In the *FTO* window alone, clone-named models
outnumber the real genes 5 to 2. The representative-transcript heuristic recovers the expected 9-exon
*FTO* structure.

**Reference-panel window extraction** (METHODS §10c). Running plink 1.9 `--r2` directly against the
67 M-variant panel makes it size its workspace from the whole file and request 34 GB, which fails; the
window is cut out with plink2 first.

## 5. Reading this study's results

There is **no independent replication cohort and no way to obtain one**. The three sample sets are
nested, so agreement between them is expected and carries little independent information. A peak here
is a threshold crossing in one scan of one sample set, and the correct next step is the GRM
random-effect follow-up in §2, not further interpretation of these *P*-values.
