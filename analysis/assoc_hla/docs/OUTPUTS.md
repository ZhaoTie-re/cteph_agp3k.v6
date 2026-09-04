# assoc_hla — output tree and data dictionary

## Read `sample_id` and every HLA id as a string

Two traps inherited from upstream, both silent:

- **107 of the 3,569 sample ids carry leading zeros** (`0000063134`, not `63134`).
  `pandas.read_csv` without `dtype` coerces them to integers and drops the zeros; an inner
  merge then loses those 107 samples with no error.
- **Marker ids contain `*` and `:`** — `HLA_A*24:02`, `AA_DQB1_57_D`. They are safe as PLINK
  variant ids and as TSV values, but not as R formula terms or shell globs without quoting.

## The tree

```
results/
├── 00.prep/
│   └── hla_gene_coords.tsv         GRCh38 coordinates, from Ensembl  [GENE_COORDS]
├── <cohort>/                       narrow_mainland | intermediate_mainland | full_mainland
│   ├── 00.prep/
│   │   ├── marker_map.tsv          every marker: id, class, gene, position, counts   [BUILD_MARKERS]
│   │   ├── marker_qc.tsv           + which were tested and why the rest were not     [MARKER_QC]
│   │   ├── pheno_covar.tsv         phenotype 0/1 and every covariate                 [PREP_PHENO_COV]
│   │   └── grm.prune.in            the genome-wide markers the GRM is built from     [PRUNE_MARKERS]
│   ├── 01.assoc/
│   │   ├── random/                 null.rda · null.varianceRatio.txt · null.fit.log  [FIT_NULL]
│   │   └── <model>/<class>/        model = fixed | random, class = allele | residue
│   │       ├── sumstats.tsv        the canonical 16-column schema + gene/position    [TO_SUMSTATS]
│   │       └── engine_extra.tsv    engine columns the schema has no slot for
│   ├── 02.omnibus/
│   │   └── omnibus.tsv             the m-1 df test at every eligible position        [OMNIBUS]
│   └── 03.signals/
│       └── significance.<model>.<class>.tsv   three thresholds, computed             [SIGNIFICANCE]
└── _run_info/                      run_manifest.json · trace · report · timeline · dag
```

## `marker_map.tsv` and `marker_qc.tsv`

One row per marker. `marker_qc.tsv` is `marker_map.tsv` plus the filter verdicts, so read the
QC file unless you specifically want the pre-filter set.

| column | meaning |
|---|---|
| `id` | the PLINK variant id — `HLA_A*24:02` or `AA_DQB1_57_D` |
| `key` | the original column name in `hla.typing`'s dosage matrix |
| `marker_class` | `allele` or `residue` |
| `gene`, `chrom`, `pos` | `pos` is the gene start plus the marker's index within the gene — a **plotting convention, not a genomic claim** (METHODS §3) |
| `position` | for residues, `<GENE>:<IMGT position>`; blank for alleles. Positions can be negative: the leader peptide runs −1 downward and there is no position 0 |
| `ac`, `an`, `af` | count of `P`, total determined chromosomes, and the ratio |
| `mac`, `maf` | the minor of `ac` / `an − ac` |
| `call_rate` | fraction of the cohort with any call at that gene |
| `determined_rate` | residues only: fraction of samples whose dosages at that **position** sum to 2. Below 1 where IMGT writes `*` or `.` — see METHODS §4 |
| `p_group` | the IPD-IMGT P group, alleles only |
| `in_reference` | **1/0 for alleles, blank for residues.** Does jMorp's 61,424 Japanese carry that P group? |
| `freq_reference` | its frequency there; blank when `in_reference = 0` |
| `frac_unconfirmed_allele` | **residues only.** Share of the marker's MINOR-state chromosomes carried by an allele the reference panel has never observed. Median 0.02; a value near 1 means the signal is an allele artefact in residue clothing — see METHODS §4 |
| `tested`, `drop_reason` | whether it entered the scan, and `mac` / `maf` / `call_rate` if not |
| `omnibus_eligible` | residue at a position with ≥ 2 tested residues and adequate determination |
| `n_residue_tested` | how many residues at that position survived |
| `arm` | `primary`, or the sensitivity arm's label |

**Two artefact flags, one per marker class, and neither is decorative.** On the first run of the
narrow cohort the single strongest fixed-model allele signal was `HLA_F*01` — a one-field call
HLA-HD could not resolve further, carried by 28 chromosomes, that a 61,424-person Japanese panel
has never observed. Every other top hit was `in_reference = 1`.

The residue equivalent earned its place the same way. After the undetermined-residue fix the
strongest residue signal was `AA_A_178_T` at *P* = 8 × 10⁻⁹ — and 15 of the 20 chromosomes
driving it carried `A*02:783`, a single allele with `in_reference = 0`. Its
`frac_unconfirmed_allele` is **0.85**.

Check `in_reference` before believing an allele result, and `frac_unconfirmed_allele` before
believing a residue one.

## `sumstats.tsv`

The canonical 16-column plink2 `--glm` schema — identical to `assoc_plink2` and `assoc_saige`,
so the same readers work — plus four columns of context.

```
#CHROM POS ID REF ALT A1 A1_FREQ TEST OBS_CT OR LOG(OR)_SE L95 U95 Z_STAT P ERRCODE
gene position cohort model
```

- **`A1` is always `P`**, and the pipeline aborts if it is not. `REF = A` means the allele or
  residue is absent, `ALT = P` that it is present, so **`OR > 1` means carrying it raises risk**.
- `L95` / `U95` are two numeric columns, never one packed string.
- `OBS_CT` is a non-missing **sample** count, the `--glm` convention, not an allele count.
- `ERRCODE` `'.'` is a usable fit. Filter on it before anything else: a Firth fit at the
  separation boundary is reported, not hidden.

## `omnibus.tsv`

One row per amino-acid position.

| column | meaning |
|---|---|
| `position`, `gene` | `DQB1:57`, `DQB1` |
| `n_residues`, `df` | residues at the position, and `df = n_residues − 1` |
| `reference_residue` | the one dropped as the reference level — the commonest, by convention |
| `residues` | all of them, `;`-separated |
| `n` | complete cases the fit used |
| `determined_rate` | see above; positions below `MinDetermined` are not here at all |
| `deviance` | 2 (log *L*<sub>full</sub> − log *L*<sub>null</sub>) |
| `P` | χ² with `df` degrees of freedom |
| `status` | `ok`, or `not_converged` / `bad_deviance` / `too_few` / `monomorphic_in_complete_cases` / `error:*` |

**Rows with `status != 'ok'` have no `P` and are not failures to hide.** A position where the
logistic fit hits the separation boundary is information about that position.

## `significance.<model>.<class>.tsv`

Three rows, one per threshold, with `primary = 1` on the one to quote.

| column | meaning |
|---|---|
| `threshold` | `effective_tests` (primary) · `bonferroni` · `genome_wide` |
| `value` | the P below which a marker is called significant |
| `basis` | how it was derived, in words |
| `n_significant` | how many tested markers clear it |
| `meff` | Li & Ji effective number of independent tests |
| `n_tested`, `n_markers_constant` | the denominator, and markers that do not vary |
| `lambda_gc` | genomic-control inflation for this scan |

Why three and not one is METHODS §6. In short: 5 × 10⁻⁸ is calibrated for a different
experiment, plain Bonferroni pretends ~1,100 near-collinear markers are independent, and the
effective count is the middle ground that is actually estimated from these data.

## What is *not* here

- **No SuSiE, no LD panels, no regional plots.** Those exist in `assoc_plink2` / `assoc_saige`
  for genome-wide scans. Fine mapping inside the MHC is what the omnibus and the conditional
  rounds do; a credible set over 150 near-collinear allele markers would not mean the same
  thing.
- **No DRB3/DRB4/DRB5 in the primary tree.** METHODS §4 gives the three reasons. Run
  `--DropGenes ''` for the sensitivity arm and read it as such.
