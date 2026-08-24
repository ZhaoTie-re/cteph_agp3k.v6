#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose : Dump the two annotation tracks a regional plot needs, for one window:
#             exons.tsv   exon/intron structure of one representative transcript
#                         per informative gene (Ensembl 86, GRCh38)
#             recomb.tsv  recombination rate from the configured bigWig track
#           Both are rendered by matplotlib downstream, so R is confined to
#           reading annotation and never to drawing — one figure style for the
#           whole component.
# Runtime : conda env `r_work`
# Used by : assoc_plink2.nf  process LD_SOURCES
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--chrom",     type = "character"),
  make_option("--start",     type = "double"),
  make_option("--end",       type = "double"),
  make_option("--recomb-bw", type = "character", default = ""),
  make_option("--out-exons", type = "character"),
  make_option("--out-recomb", type = "character")
)))

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "gene_utils.R"))

chr_plain <- sub("^chr", "", opt$chrom)

# ── gene models with exon structure ────────────────────────────────────────
exons <- representative_exons(informative_genes(opt$chrom, opt$start, opt$end))
write.table(exons, opt$`out-exons`, sep = "\t", quote = FALSE, row.names = FALSE)

# ── recombination rate ─────────────────────────────────────────────────────
rec_df <- data.frame(start = numeric(), end = numeric(), rate = numeric())
if (nzchar(opt$`recomb-bw`) && file.exists(opt$`recomb-bw`)) {
  try({
    suppressPackageStartupMessages(library(rtracklayer))
    # The peak table carries plink's #CHROM ('16'); this bigWig is named 'chr16'.
    # A wrong seqname returns silently empty, which would look like a region with
    # no recombination data, so try every spelling.
    avail <- names(seqinfo(BigWigFile(opt$`recomb-bw`)))
    cand <- c(opt$chrom, chr_plain, paste0("chr", chr_plain))
    hit <- cand[cand %in% avail]
    if (length(hit)) {
      gr <- import(opt$`recomb-bw`,
                   which = GRanges(hit[1], IRanges(opt$start, opt$end)))
      if (length(gr) > 0)
        rec_df <- data.frame(start = start(gr), end = end(gr), rate = gr$score)
    }
  }, silent = FALSE)
}
write.table(rec_df, opt$`out-recomb`, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("[region_tracks] %s:%.0f-%.0f  genes=%d  exons=%d  recomb_intervals=%d\n",
            opt$chrom, opt$start, opt$end,
            length(unique(exons$gene)), nrow(exons), nrow(rec_df)))
