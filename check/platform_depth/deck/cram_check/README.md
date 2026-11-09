# cram_check: reads at the two flagged leads (deck, 2026-11-06)

User request: "在27和28，我需要你查看cram原文件来判明问题（各自新加页面）". Read only; every full-mainland genome
(3,101 CRAMs, the indexes registered in results/00.samples/sample_table.tsv; reference hs38DH.fa).

- `cram_scan.py`: one row per genome.
  - ADARB2 (chr10:1747167 CGTG>C): the lead lies in a minisatellite of about 390 bp, 31-bp units. The carrier haplotype is
    one unit shorter, called as a 3-bp and a 28-bp deletion. The aligner places its reads without a gap (one unit
    shifted, a few mismatches), so the CIGAR does not show it. Every read of the minisatellite (any MAPQ; unmapped
    mates placed there included) is searched for the carrier signature GATCGCATGATGTGTA, which is absent from the
    reference window, and for its reference form GATCGCGTGATGTGTA (three units).
  - LINC01146 (chr14:88019710 G>A): G and A bases, by strand, with MAPQ and base quality >= 20.
- `signature_discovery.py`: how the signature was found. It ranks k-mers by how well they separate HiSeq X controls
  called 1/1 (8) and 0/1 (25) from 0/0 (40), using reads dumped by `dump_reads.py` (any MAPQ).
  It ranks K = 15, 19 and 21; the 16-mer used, GATCGCATGATGTGTA, spans the junction that the two called deletions
  create (the reference window with both deletions applied holds it; no reference unit does), and its reference form
  GATCGCGTGATGTGTA sits in three units (chr10:1747131, 1747162, 1747193). The minisatellite spans chr10:1746881–1747271
  (`stats.json` `cram.vntr`), 12–13 units of 31 bp.
- `merge_calls.py`: the join to the sample table and to the pre-QC calls at the two leads
  (wgs.auto.par/results/03_annotate_af_norm_gt, bcftools query).
- `cram_scan.tsv`: the per-genome table read by collect_stats.py (key `cram`) and make_figures.py
  (cram_ADARB2, cram_LINC01146). Internal: holds sample IDs; slides show counts only.

A read-level carrier at ADARB2 holds at least 4 signature reads. The counts are bimodal: 0–2 reads are noise; called
carriers have 4–20 reads (185 of 186 controls). With 3 instead of 4, no T7 count changes.
- `recall_check.py` (2026-11-06, the restructured flag slides): where the ADARB2 deletion is lost. Writes `aligner.tsv`
  (every full-mainland CRAM's alignment record, the bwa command normalized: program and reference by base name, thread
  count and read group dropped) and `recall_adarb2.tsv` (the genotype at the lead after Q1 re-called the deep cases'
  reads at full depth, baseline arm of results/03.calls/downsampling/arms.qc.vcf.gz, joined to cram_scan.tsv).
  Read only; internal (sample IDs). Usage: `python3 -B recall_check.py [n_proc]`.
- `hc_debug/` (2026-11-06, user: "我们先看看HaplotypeCaller 的诊断文件"): the cohort's HaplotypeCaller command re-run on
  14 genomes at the ADARB2 lead with output-only diagnostics (`run_hc_debug.sh`), the assembly log summarized
  (`parse_hc_log.py` → `hc_summary.tsv`) and the tables the deck reads (`hc_tables.py` → `hc_genotypes.tsv`,
  `cohort_records.tsv`, `cohort_ga.tsv`; stats key `cram.hc`). See `hc_debug/README.md`. Internal (sample IDs).
- `hc_census/` (2026-11-07, user: "我们有没有什么办法取样取的更加系统有说服力，另外把样本ID也可以表上"): the same command on
  every full-mainland carrier in the reads or the calls (292) and 20 random non-carriers per platform (392 genomes; no
  selection). The deletion is assembled exactly where it was called (392 of 392); it replaces the hand-picked pilot in
  the deck (panel b of slide 26; stats key `cram.hc_census`, counts only). See `hc_census/README.md`. Its per-sample
  table holds sample IDs and is git-ignored.
- `tommo_signature.py` → `tommo_signature.tsv` (2026-11-07, round 12 item 10; allele counts only): the ToMMo 61KJPN
  frequencies of the two alleles that give the ADARB2 signature, the one-unit deletion (AF 4.73%) and the G>A in a
  full-length unit (1747168: 0.036%; 1747199: 0.043%, VQSR-filtered), with the carrier share each predicts under
  Hardy–Weinberg. Stats key `cram.tommo_sig`.
