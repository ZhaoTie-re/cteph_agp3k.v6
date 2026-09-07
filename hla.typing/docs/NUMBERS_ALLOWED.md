# Numbers a document in this component may state without deriving them

This file is read by `verify.sh` §9. It is a `.md` and not a `.txt` because the
repository's `.gitignore` is an extension allow-list with no `.txt` rule, so a
`.txt` version of this file never reached the remote and a fresh clone could not
run its own checks.

```
# Numbers a document in this component may state without deriving them.
#
# verify.sh section 9 fails on any number in a .md that is neither in
# results/_run_info/facts.json nor listed here. Everything here is a CONSTANT --
# a version, a parameter, a fixed property of an external resource, a section
# number -- not a result. If a value moves when the data move, it does not
# belong in this file; it belongs in build_facts.py.
#
# One entry per line. Text after # is a comment.

# ---- pinned software and reference releases ----
1.7.1        # HLA-HD
3.64.0       # IPD-IMGT/HLA dictionary and protein alignments
3.50.0       # HLA_gene.split, pinned beside the dictionary
2.5.5        # bowtie2
1.0          # nextflow manifest version

# ---- parameters, which are chosen and not measured ----
0.95         # hlahd.sh -c, the cutting rate
100          # hlahd.sh -m, the minimum read length
2            # AlleleFieldDepth, and MinNumVar-style small counts
1            # MinResidueCount
0            # counts that are asserted to be zero by design
20           # MIN_CHR in plot_allele_freq.py
0.98         # R_GOOD, the concordance a locus must clear
0.02         # MIN_DIFF_TO_NAME
17           # the matched depth window, lower bound
21           # the matched depth window, upper bound
7.20         # figure width in inches
600          # figure dpi
4571         # MB per core on this site's scheduler

# ---- fixed properties of the MHC window and the references ----
6            # chr6
29540000     # MhcStart
33420000     # MhcEnd
33           # genes HLA-HD types
19           # genes IMGT publishes a protein alignment for
13           # loci the primary reference carries
5            # loci the secondary reference carries
150          # the N150 dictionary suffix

# ---- typography and section numbering ----
1            # section numbers, list numbers, field depths
2
3
4
5
6
7
8
9
10
95           # 95 % confidence intervals
2026         # the year of this work
2023         # Tadaka et al., the jMorp citation
2018         # the 1000 Genomes HLA panel release
2025         # dates in provenance notes

# ---- fixed properties of the inputs and the site, not results ----
543          # extraction regions in contig_list.txt (chr6 window + HLA decoys)
111          # CRAMs in this study with no index anywhere, built by INDEX_CRAM
88           # TB of CRAM read by a full run
480          # CPU-hours for a full run
46,005       # alleles in the pinned IPD-IMGT/HLA 3.64.0 dictionary
3,366        # sequences in hs38DH.fa, matching the CRAM @SQ
525          # HLA* decoy contigs in hs38DH
67           # the largest CRAM, in GB
48           # maxForksExtract, the concurrency cap on CRAM reading
1.2          # GB of HLA-HD output per sample
90           # MB of extracted FASTQ per sample
98           # per cent of HLA-HD output that is intermediate SAM
17           # PopGMM components forming the mainland cohort
36,568       # MB requested by HLAHD
73,136       # MB on the HLAHD retry
18,284       # MB for the small processes
8,000        # MB for INDEX_CRAM

# ---- the same constants as written in prose, with separators ----
29,540,000   # MhcStart
33,420,000   # MhcEnd
9,142        # MB for INDEX_CRAM after the 4571 correction
18,284       # MB for the small processes
36,568       # MB requested by HLAHD
73,136       # MB on the HLAHD retry
3.6          # GB, the low end of HLAHD's measured peak RSS
6.5          # GB, the high end
5.3          # GB
1.51         # GB, peak RSS of one measured sample
1.17         # GB of alignments per sample
4.2          # TB of alignments across the cohort
321          # GB of extracted reads across the cohort
4.5          # TB, the cost of copying instead of symlinking
34           # MB, the residue reference
100          # MB, everything that is copied
200          # seconds saved by not re-reading the FASTA per task
30           # minutes, the reference implementation's per-sample time
50           # minutes, the same, upper bound
239,123      # symlinks the naive COLLECT_ALLESES staging would have made
16           # chr6_*_alt contigs
99.622       # per cent of T7 reads surviving -m 100
99.628       # the same for HiSeqX
12.0         # x, the shallowest measured MHC coverage
0.3          # per cent, the rarest allele a small panel can sample
580          # the size ratio between the two reference panels
0.034        # the standard error near a frequency of 0.4 on the small panel

# ===========================================================================
# ONE-OFF MEASUREMENTS, recorded once and never recomputed.
#
# Everything above this line is a constant. Everything BELOW it is a number
# that was measured once — during design, validation or a debugging session —
# and written down because the reasoning depends on it. They are not derived
# from results/ and a re-run does not update them, which is exactly why they
# have to be listed here rather than silently trusted: a reader who wants to
# know whether a number is live checks whether it is in facts.json, and if it
# is in this block instead, the answer is no.
#
# If one of these ever needs to be current, the fix is to derive it in
# build_facts.py and delete the line here.
# ===========================================================================

# --- IMGT alignment structure and the P-group translation (METHODS 11, 13) ---
25,290       # two-field names in hla_nom_p.txt, asserted single-valued
326,922      # the global registry count behind HLA-HD's freq_data prior
120,000      # rows a per-(sample, allele) long form would produce
151          # exon-groups affected by IMGT's exon-boundary convention
461          # P groups containing more than one 2-field allele
460          # of those whose members differ only outside the ARD
99.8         # per cent, the same as a share
22,572       # residue differences between members of a P group
21,717       # of those outside the antigen recognition domain
855          # of those inside it
96.2         # per cent outside the ARD
21,700       # residue differences P-group collapsing would delete
76           # per cent of DRB4 chromosomes carrying DRB4*01:03
0.5430       # frequency of the DRB4 null allele in the reference
0.625        # the low end of r before the P-group translation
0.991        # the high end
186          # chromosomes the 1000 Genomes panel types at DQB1, of 210

# --- HLA-HD behaviour, measured during design (METHODS 4, 8) ---
348          # samples in the reference cohort used for these measurements
359          # heterozygous calls in one measured implementation
6            # samples that failed silently, of 348
232          # samples in the reference cohort for the DRB3/4/5 measurement
280          # of 348 samples showing the DRB3/4/5 pattern
13.5         # per cent of HLA-A calls showing '-' where it means homozygous
117          # DRB4 'one allele + -' calls in the reference cohort
163          # DRB5 'Not typed' calls in the reference cohort
61           # DRB3 'one allele + -' calls in the reference cohort
21           # a small count in the same table
3            # the same
45           # bases of a 150 bp read matching, against 30 of a 100 bp read

# --- the two panels' agreement where they overlap (OPEN_QUESTIONS 2) ---
0.986        # per-locus r between the panels, low end
0.988
0.998
0.4018       # our DRB3 non-null frequency
0.4565       # the reference's
0.674        # reference DRB3*01:01P frequency
0.180        # ours

# --- the deferred accuracy experiment (OPEN_QUESTIONS 2) ---
2,693        # rows in the 1000 Genomes HLA type file
104          # JPT samples with both a published type and a downloadable CRAM
115          # CRAMs in scope for the experiment
80           # 1KG CRAMs already on this system
0.39         # the share of them that are JPT
15.1         # GB per CRAM
1.57         # TB to download
10.4         # hours at the measured rate

# --- residue distance between a tail call and its family (OPEN_QUESTIONS 3) ---
# Measured on the superseded 3,569-sample run; a property of the IMGT
# alignments, not of the cohort. Carried over deliberately, see that section.
396          # A: tail-allele pairs compared
390          # A: positions in the alignment
338          # C: pairs
417          # C: positions
395          # B: pairs
483          # B: positions
117          # DQB1: pairs
47           # DQB1: median residues differing
274          # DQB1: positions
144          # DRB1: pairs
85           # DRB1: median residues differing
290          # DRB1: positions
51           # per cent of A pairs differing by <= 2
55           # per cent of C pairs
41           # per cent of B pairs
15           # per cent of DQB1 pairs
22           # per cent of DRB1 pairs

# --- residue matrix structure measured once (OUTPUTS) ---
2,703        # alignment rows where HLA-A position -22 is '*'
201          # alleles where the counts appeared to disagree
1,543        # alleles in that comparison
```
