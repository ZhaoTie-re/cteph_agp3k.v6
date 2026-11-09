# hla.typing.v2 — HLA genotypes and the evidence for their quality

**What it is.** HLA genotypes of the full_mainland cohort (3,101 samples: 2,662 AGP3K controls, 439 PH cases)
at 19 HLA genes, each call with its read evidence and a confidence class, and a quality summary per locus
against two Japanese reference panels.

**What it is not.**
- **Not a new typing.** HLA-HD 1.7.1 ran once over these samples in `../hla.typing` (2026-08-25); its output
  is archived here and everything else is derived from it.
- **Not a study of platform or depth.** Those are left to a later down-sampling and platform-sensitivity
  analysis.
- **Not an association study.**

Methods, rules and thresholds: [docs/METHODS.md](docs/METHODS.md). Files: [docs/OUTPUTS.md](docs/OUTPUTS.md).
Figures: [results/figures/README.md](results/figures/README.md). Slides: [deck/](deck/).

## Run

```
./run.sh              # all steps, 00a → 06 (minutes, except 00a: ~40 min of Lustre reads, once)
./run.sh 03 04 06     # selected steps
./verify.sh           # checks every output against the rules in docs/METHODS.md
```

Plain Python scripts, run in order; there is no workflow engine, because nothing here is expensive.
Interpreter: `cteph_geno_pro` (absolute path in `run.sh`). Paths, versions and every threshold are in
`scripts/config.py`.

| step | script | does |
|---|---|---|
| 00a | `00_archive_inputs.py` | one `tar.gz` per sample: HLA-HD results, run log, validation, task scripts |
| 00b | `00_provenance.py` | checks one version, one command and validation `ok` for all samples; read counts; md5 of every pinned file |
| 01 | `01_calls.py` | HLA-HD candidates → 2-field genotypes, states, *DRB3/4/5* copy number from *DRB1* |
| 02 | `02_evidence.py` | reads and unique reads per allele, flags, call confidence, sample QC |
| 03 | `03_panels.py` | controls vs ToMMo 61KJPN and 1KG JPT; homozygosity and dropout |
| 04 | `04_confidence.py` | locus tiers |
| 05 | `05_exports.py` | allele and residue dosage tables (19 genes, every sample); a residue the IPD-IMGT/HLA alignment leaves undetermined is blank |
| 06 | `06_assoc_inputs.py` | **the matrices the association reads**, one set per cohort: 16 genes, call rate ≥ 95 % in cases and in controls, then MAF ≥ 0.5 %; plus every column's measures and verdict. No test is run |
| 07 | `07_figures.py` | eight figures |

## Layout

```
hla.typing.v2/
  run.sh  verify.sh
  scripts/          config.py · hlahd_io.py · nomenclature.py · panels.py · figlib.py · 00–07 steps · build_residue_ref.py
  assets/           fonts/ (Liberation Sans, OFL) · imgt/Allelelist_history.txt (IPD-IMGT 3.64.0) ·
                    pdb/ (1HHK, 1DLH) and ucsc/ (chr6 bands), for the deck
  results/
    00.inputs/      typing/<sample>.tar.gz · archive_manifest.tsv · sample_inputs.tsv · provenance.json · hlahd_freq_data/
    01.calls/       calls_long.tsv · genotypes.tsv · tie_candidates.tsv · drb_reconciliation.tsv
    02.evidence/    call_evidence.tsv · gene_summary.tsv · sample_qc.tsv · drb_consistency.tsv
    03.panels/      locus_concordance · allele_concordance · homozygosity · panel_ceiling · shift_direction ·
                    unconfirmed_alleles · confidence_vs_panel · excluded (.tsv)
    04.confidence/  locus_tiers.tsv
    05.exports/     allele_dosage · residue_dosage · allele_pgroup_map · call_confidence (.tsv) · residue_reference/
    06.assoc_inputs/ <cohort>/allele_dosage.tsv · <cohort>/residue_dosage.tsv · markers.tsv · qc_record.tsv · summary.json
    figures/        01–08 .png · README.md
  docs/             METHODS.md · OUTPUTS.md
  deck/             the slides: hla_typing.pptx (24 slides) and their sources; every marker number is read from 06.assoc_inputs
```

## Results: the basic evidence that the genotypes are usable

Two checks, one figure (`results/figures/03.panel_frequencies.png`):

1. **Complete.** All 3,101 genomes typed at 19 HLA genes and every run validated. Of 58,919 calls,
   99.9 % have a result and 99.0 % a 2-field genotype; the lowest gene is *DRB4*, at 91.5 %.
2. **Consistent with the Japanese population.** Allele frequencies in the 2,662 controls agree with those of
   61,424 Japanese individuals in ToMMo 61KJPN at all 11 compared loci (*A, B, C, E, F, G, DRB1, DQA1, DQB1,
   DPA1, DPB1*), with Pearson r = 0.975–1.000, each r computed within one locus. Not compared, and why:
   - *DRA, DMA, DMB, DOA, DOB, DRB5*: ToMMo publishes no frequencies for them (nor does 1KG JPT);
   - *DRB3*: ToMMo assigns a *DRB3* allele to every chromosome, though only DR52 haplotypes carry it;
   - *DRB4*: its presence follows *DRB1*, which is already compared; ToMMo's *DRB4* entries do not match the
     DR53 share of its own *DRB1*; and the allele mix among carriers is dominated by *DRB4\*01:01*.
   *DRB3/4/5* are also left out of the association analysis (see the deck).

For a manuscript:

> HLA alleles were typed from the whole-genome sequences with HLA-HD (Kawaguchi et al., 2017). Allele
> frequencies in the 2,662 controls were highly concordant with those of 61,424 Japanese individuals in the
> ToMMo 61KJPN panel (Pearson r = 0.975–1.000 across 11 loci).

**Supplementary, if asked.** The further checks stay in `results/` and `results/figures/`; they are not in the
deck and not needed for the basic claim:
- the reads behind each call, and the confidence classes;
- the *DRB3/4/5* copy-number reconciliation;
- the comparison with 1KG JPT;
- the locus tiers.

`analysis/assoc_hla` has not been re-run on these exports. The *DRB3/4/5* dosages changed: a single copy is
now 1 and an absent gene 0. Under the association's marker rule (a residue left undetermined by IPD-IMGT/HLA set
missing; minor allele count ≥ 20; call rate ≥ 95 %), the 16 obligate genes give 163 alleles and 697 residues in
the full mainland cohort (146–163 and 682–697 across the three cohorts, every genome typed); the deck's last two
slides list these inputs and the recommendations for the association analysis.
