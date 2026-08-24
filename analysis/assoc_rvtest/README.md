# assoc_rvtest

Gene-based rare-variant association (RVTESTS) over three nested cohorts, on the call set and
the minor-allele-count threshold that [`tuning.rv`](../../tuning.rv/README.md) has already
decided.

→ Method and rationale: **[docs/METHODS.md](docs/METHODS.md)** ·
Output tree and data dictionary: **[docs/OUTPUTS.md](docs/OUTPUTS.md)** ·
This study's numbers: **[docs/STUDY_NOTES.md](docs/STUDY_NOTES.md)** ·
Shared visual grammar: **[../_shared/docs/FIGURES.md](../_shared/docs/FIGURES.md)**

## What this component does *not* establish

The design is **fully confounded**: cases are 30× (DNBSeq/NovaSeq), controls 15× (HiSeqX),
with zero platform overlap. `tuning.rv` quantifies that and returns the minAC at which the
*adjusted* apparent group effect reaches a statistical null. Inheriting its threshold is what
makes a burden test interpretable here at all. **A gene-based hit is a candidate, not a
confirmed association**, and there is no independent cohort for this phenotype.

## The contract with `tuning.rv`

This component derives none of its QC. Per cohort it reads:

| what | from | value in this study |
| --- | --- | --- |
| sample- and variant-QC'd rare call set | `02.callset_filter/filtered.{bed,bim,fam}` | 14.75 M / 16.51 M / 18.95 M variants |
| the minAC to apply | `05.qc_collect/minac_recommendation.tsv` | **2** in all three, `calibrated=true` |

`params.MinAC = null` means *inherit*. Set it to a number only to answer a "what if"
question — the run manifest records which of the two happened.

If `tuning.rv` has not run for a cohort, this pipeline **fails at launch with the path it
wanted**, rather than silently analysing an unfiltered call set.

## Run it

```bash
source activate dsl2
cd analysis/assoc_rvtest
nextflow run assoc_rvtest.nf -resume
nextflow run assoc_rvtest.nf --Cohorts narrow_mainland          # one cohort
nextflow run assoc_rvtest.nf --MinAC 5                          # override the tuning value
```

Analysis scripts run under the `cteph_geno_pro` conda env, activated inside each process.

### Parallelism

Everything from `SPLIT_CALLSET` to `RVTEST` runs on **balanced bins of whole chromosomes**.
No gene spans a chromosome and rvtest tests genes independently, so the split is **exact** —
the concatenated result is identical to a whole-genome run, not an approximation.

The one thing that cannot be split is the Benjamini-Hochberg correction, which is a property
of the whole scan. `MERGE_ASSOC` is the barrier that guarantees it: the per-bin tables become
one table before any tier is assigned.

| | before | after |
| --- | --- | --- |
| critical path, per cohort | 203 min | **~19 min** |
| RVTEST jobs | 18 | 216 |

Twelve bins, not twenty-two: `chr2` alone carries 9.05 % of variants and is the floor no
split can beat, and twelve bins already reach it (max bin 9.05 %) while running 45 % fewer
jobs than one-bin-per-chromosome. `params.ChromBins` holds the packing.

### Key parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `Cohorts` | all three | cohorts to analyse (comma list on the CLI) |
| `MinAC` | `null` | `null` = inherit from `tuning.rv`; a number overrides |
| `ChromBins` | 12 balanced bins | the unit of parallelism |
| `impactFilters` | L+M+H, M+H, H | snpEff IMPACT strata |
| `testMethods` | `skato`, `cmc` | one kernel test and one burden test |
| `MinNumVar` | 3 | minimum qualifying variants for a gene to be tested |
| `FdrSignificant` / `FdrSuggestive` | 0.05 / 0.10 | the two tiers |
| `PcLabel` / `NPcs` | `bbj_mainland` / 10 | ancestry PCs in the covariate set |
| `PeakStratum` / `PeakMethod` | `moderate_high` / `skato` | the scan the cross-cohort figure draws |

## Tiers

`significant` = BH FDR < 0.05, `suggestive` = BH FDR < 0.10, computed **per scan**
(one scan = one cohort × stratum × method).

FDR rather than a fixed *P* because a scan tests 164 genes under HIGH and ~8.6 k under
MODERATE+HIGH; a fixed cut would mean a different multiple-testing burden per stratum. The
price is that **the implied *P* differs between strata** — measured here, FDR 0.05 falls at
1.3 × 10⁻⁴ in HIGH and 4.5 × 10⁻⁵ in MODERATE+HIGH. Every scan figure therefore draws *its
own* BH-implied lines and prints their values; a tier label is comparable within a stratum,
not across strata.

## Scripts

```
scripts/
  rvtest_prepare_{main,tools}.py   PLINK -> VCF + rvtest-format pheno/covar/refFlat
  snpeff_anno_{main,tools}.py      snpEff annotation
  info_filter_{main,tools}.py      subset a VCF by predicted IMPACT
  rvtest_post_process.py           minimum-NumVar filter
  gene_scan.py                     genome-wide BH + tier assignment + scan QC
  plot_gene_scan.py                gene-level Manhattan + QQ
  gene_detail.py                   per-gene deep dive (variants, carriers, ToMMo)
  plot_gene_detail.py              one figure per tiered gene
  cross_cohort_genes.py            the union of tiered genes, in every cohort
  plot_cohort_compare_genes.py     the three cohorts side by side
../_shared/scripts/
  plot_style.py                    style, tier lines, layout helpers
  figure_doc.py                    a standalone figure's .md, a fan-out figure's stats.json
  figure_catalogue.py              one document per figure family per cohort
```

## Checks worth running after a run

```bash
# the minAC that was actually applied, and where it came from
python3 -m json.tool results/_run_info/run_manifest.json | grep -A6 minac

# nothing failed
awk -F'\t' 'NR>1 && $4!="COMPLETED" && $4!="CACHED"' results/_run_info/trace.txt

# both tiers populated, and the figure count matches the gene count
for c in narrow_mainland intermediate_mainland full_mainland; do
  for t in significant suggestive; do
    printf '%s %s: %s genes, %s figures\n' "$c" "$t" \
      "$(ls results/$c/05.genes/$t 2>/dev/null | wc -l)" \
      "$(ls results/$c/figures/02.gene/$t/*.png 2>/dev/null | wc -l)"
  done
done

# every fan-out family has its catalogue; standalone figures have their sidecar
find results -name '*.png' | while read p; do
  d=$(dirname "$p")
  case "$d" in */significant|*/suggestive) [ -f "$d/../README.md" ] || echo "NO CATALOGUE $p" ;;
                *) [ -f "${p%.png}.md" ] || echo "MISSING $p" ;; esac
done
```

The covariate assertion is not one of these: it runs **inside** `RVTEST`, because rvtest
exits 0 when a `--covar-name` column is missing and silently fits an uncovaried model. A run
that reaches the end has had its covariates checked 216 times.
