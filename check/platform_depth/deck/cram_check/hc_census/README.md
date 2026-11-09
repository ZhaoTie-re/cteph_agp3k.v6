# hc_census: HaplotypeCaller at the ADARB2 lead in every carrier (2026-11-07)

User request: "我们有没有什么办法取样取的更加系统有说服力，另外把样本ID也可以表上". The pilot (`../hc_debug`) diagnosed 14
genomes chosen by hand; the census diagnoses every carrier, so no selection is involved. Read only on the CRAMs.

**Sample IDs.** `samples.tsv`, `hc_census.tsv`, `out/` and the SLURM log hold sample IDs and are git-ignored (the
repository's extension allow-list has no rule for them). This README, the scripts and `stats.json` hold none. The deck's
appendix slide lists the case carriers by ID; the deck PDF is excluded from git for that reason (`.gitignore`).

## Genomes (`make_samples.py`)

Full mainland, which holds narrow and intermediate mainland: every genome that carries the deletion in its reads (at
least 4 signature reads, `../cram_scan.tsv`) or in the cohort's calls (0/1 or 1/1 before QC), and, as negative controls,
20 genomes per platform drawn at random (seed 20261107) from the clear non-carriers (no signature read, called 0/0).

| Group | Genomes |
|---|---|
| Called carriers, 150-bp reads | 190 (185 HiSeq X controls; 5 cases: NovaSeq 2, G400RS 15x 2, G400RS 30x 1) |
| Called carrier without signature reads | 1 (HiSeq X) |
| Read carriers the calls missed, 150-bp reads | 69 (65 HiSeq X controls; 4 cases: G400RS 15x 3, NovaSeq 1) |
| DNBSeq-T7 read carriers (100-bp reads) | 32 (all missed by the calls) |
| Random non-carriers | 100 (20 per platform) |
| Total | 392 |

## Run (`run_hc_census.sh`)

The pilot's command unchanged (GATK 4.1.4.0, Java 8, the cohort's `CALL_GVCF` settings, `-ERC GVCF`, calling interval
chr10:1735001-1760000, reads of the extraction window chr10:1734001-1761000), with `--debug-assembly` and a bamout of
every assembled haplotype (`ALL_POSSIBLE_HAPLOTYPES`). SLURM job 26577887, queue eb, `--rsc p=1:t=16:c=16:m=73136M`,
8 genomes in parallel, 12 min 10 s; no error in any log.

## Tables (`census_tables.py`)

- `hc_census.tsv` (internal: sample IDs): per genome, the assembly region over the lead (reads, k-mer sizes, candidate
  haplotypes); whether a candidate haplotype removes exactly one 31-bp unit inside the minisatellite (deleted minus
  inserted bases within chr10:1746881-1747271 equal 31, wherever the gap is placed); the genotypes the re-run writes at
  the lead (1747167), the 28-bp deletion (1747171) and the signature's G>A positions (1747168, 1747199); whether the
  re-run reproduces the cohort's call at the lead (any ALT index counts as a carrier, so 1/0 does too).
- `hc_census_groups.tsv`: the same by group, counts only.

## Findings

1. The re-run reproduces the cohort's call at the lead in 392 of 392 genomes.
2. The deletion haplotype is assembled exactly where the deletion was called (392 of 392): in all 190 called carriers
   with 150-bp reads and the one without signature reads; in none of the 69 read carriers the calls missed on 150-bp
   reads, none of the 32 T7 carriers, and none of the 100 random non-carriers.
3. Every assembly over the lead uses k = 85 (smaller k dropped for non-unique reference k-mers), as in the pilot.
4. Read support does not separate the T7 carriers: their median signature reads (8) sit between the called (9) and the
   missed (6) HiSeq X carriers, yet none is assembled; on 150-bp reads the missed carriers have fewer signature reads.
5. A heterozygous G>A at the signature's position (1747168 or 1747199) is written instead in 17 of the 32 T7 carriers
   (GQ 6-99) and 37 of the 65 missed HiSeq X carriers; in 3 called carriers; in no non-carrier. The pilot's "every T7
   carrier" held for its 5 genomes only.
6. One missed NovaSeq carrier (4 signature reads) has no assembly region over the lead: its miss is an inactive region,
   not a failed assembly.

The read-length cause stays inferred: the census shows where the deletion is lost (assembly) and that it is lost in
every T7 carrier, not why. The direct test (150-bp carrier reads cut to 100 bp and re-called) was not run (user
2026-11-07: "先不做").

## Use in the deck

`collect_stats.py` key `cram.hc_census` (counts only); `make_figures.cram_adarb2` panel b (share reconstructed per group,
95% Wilson CI); slide 26 step 4 and notes; appendix slide `apHcIds` (case carriers by sample ID, read from
`hc_census.tsv` by `build_slides.js`; `deck_check.py` allows IDs on that slide only and takes its per-sample numbers from
`hc_census.tsv`).
