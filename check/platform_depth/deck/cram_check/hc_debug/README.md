# hc_debug: HaplotypeCaller diagnostics at the ADARB2 lead (2026-11-06)

User request: "我们先看看HaplotypeCaller 的诊断文件". Read only on the CRAMs; internal (sample IDs).

**Superseded in the deck (2026-11-07) by `../hc_census`**, which runs the same command on every carrier instead of
14 genomes chosen by hand; the pilot's finding 5 ("all 5 T7 carriers" written as a G>A) holds in the census for 17 of
32 T7 carriers.

## What was run

`run_hc_debug.sh` (SLURM job 26552941, queue eb, 70 s): for the 14 genomes of `samples.tsv`, the reads of the
extraction window chr10:1734001-1761000, then the cohort's HaplotypeCaller command exactly as `CALL_GVCF` of
`platform_depth.nf` runs it (GATK 4.1.4.0, Java 8, `-ERC GVCF`, defaults, calling interval chr10:1735001-1760000 of
`results/02.targets/hc.intervals`), with output-only diagnostics: `--debug-assembly` (`<id>.hc.log`), `--bam-output`
with `ALL_POSSIBLE_HAPLOTYPES` (`<id>.bamout.bam`), `--assembly-region-out`, `--activity-profile-out`, `--graph-output`.
`parse_hc_log.py` summarizes the assembly region over the lead per genome (`hc_summary.tsv`).

Genomes: 5 DNBSeq-T7 read carriers (10–15 signature reads), 1 T7 non-carrier, 3 HiSeq X called carriers, 2 HiSeq X
read carriers the cohort missed, and one called carrier each on G400RS 30x, G400RS 15x and NovaSeq (150-bp reads;
G400RS shares the T7's DNB chemistry).

## Findings

1. The calls reproduce the cohort's: no T7 carrier gets the 3-bp + 28-bp deletion; the 150-bp called carriers do.
2. The minisatellite is active in every genome (assembly regions cover it), so the miss is not an inactive region.
3. In the assembly region over the lead, every genome drops k = 10, 25, 35, 45, 55, 65 and 75 because the reference
   window holds non-unique k-mers (the 31-bp repeat) and assembles with k = 85. A 100-bp read holds 16 such k-mers,
   a 150-bp read 66.
4. Assembled haplotypes: the deletion haplotype (CIGAR with 3D…28D) appears in all 6 called 150-bp carriers (3 HiSeq X,
   1 each G400RS 30x, G400RS 15x, NovaSeq; 1–3 haplotypes each), in none of the 5 T7 carriers and in neither of the 2
   missed HiSeq X carriers. (Corrected 2026-11-07: an earlier version of this note said "5 called 150-bp carriers".)
5. Without it, the caller writes the carrier's one distinct base as a G>A SNV inside a reference unit: chr10:1747168
   or chr10:1747199 (the G of the signature's reference form in units 2 and 3), heterozygous with GQ 99 at one of the
   two in all 5 T7 carriers (e.g. AD 8,13; hc_genotypes.tsv).
6. In the cohort's calls (wgs.auto.par 01_prepare_vcf, VQSR PASS) chr10:1747168:G:A exists (AC 39) and is carried only
   by read carriers whose deletion was not called: T7 10 of 32, HiSeq X 22 of the 65 missed, G400RS 15x 1 of 3;
   none of the 191 called carriers (and 3 non-carriers; cohort_ga.tsv). AC 39 counts the cohort VCF's 3,592 genomes; 36 of
   the carriers are in full mainland. It was removed by the cohort's site QC (02_filter_vqc: `VQSLOD > 10 & MQ > 58.75`;
   its VQSLOD 6.15; the deletion's 10.15). chr10:1747199:G:A is absent from the cohort's calls.

So the T7 carriers' haplotype is not assembled as the deletion; part of it surfaces as a G>A SNV in the repeat, and
that record fails the site QC. The k = 85 / read-length explanation is inferred (finding 3), not yet tested: candidate
tests are `--min-pruning 1` or a smaller k with `--allow-non-unique-kmers-in-ref` on the T7 carriers, and 150-bp
carrier reads cut to 100 bp.

## Tables (hc_tables.py, 2026-11-07)

`python3 -B hc_tables.py` (cteph_geno_pro; pysam) writes, read only:
- `hc_genotypes.tsv` (internal, sample IDs): per diagnosed genome, the deletion haplotypes assembled, the k-mer size,
  and the genotype written at 1747167, 1747171, 1747168 and 1747199 in its debug GVCF.
- `cohort_records.tsv`: the cohort's records in chr10:1747160-1747205 before the site QC (01_prepare_vcf) with VQSLOD,
  MQ, AC, AN, and whether each survives it (02_filter_vqc).
- `cohort_ga.tsv` (counts only): the full-mainland genomes of ../cram_scan.tsv by group, read carriers called and missed,
  and how many carry the G>A at 1747168 before the site QC.
The deck reads them through collect_stats.py (key `cram.hc`).
