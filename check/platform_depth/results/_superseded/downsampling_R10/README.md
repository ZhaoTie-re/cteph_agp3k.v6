# Superseded: the R = 10 down-sampling of the 89 deep cases (2026-10-30 to 2026-11-03)

The first full run down-sampled the 89 NovaSeq and DNBSeq-G400RS 30x cases only, each to the HiSeq X controls'
median effective depth, with a baseline and R = 10 replicates (run full.1; report.3 and report.4).

On 2026-11-03 the user set the design to every case platform but DNBseq-G400RS 15x at the controls' median depth
(the 319 DNBSeq-T7 cases added with one platform fraction, D* / their median = 0.96528) and R = 3 replicates for
every platform, for uniformity. The new run reuses the 89 deep cases' baseline and r01-r03 from the cache; r04-r10
are no longer part of the analysis.

This folder keeps the R = 10 record as it was published: the arms' genotype, HLA and re-test folders (baseline,
r01-r10), the down-sampling summaries, the arm plan and joint calls, the gates G6, G8 and G9, and the report tables
of that run (11.array, 13.followup, 14.evidence as of report.4: full-mainland arrays and follow-up). Nothing reads it.
