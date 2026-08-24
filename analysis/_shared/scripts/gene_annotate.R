#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose : Assign a gene to each lead variant for the `Gene` column of the
#           lead-variant table. Reads a TSV with chrom/pos, writes chrom/pos plus
#           Gene, Gene_Biotype and Gene_Distance_bp.
# Runtime : conda env `r_work`
# Used by : scripts/annotate_leads.py (shells out to this)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--positions", type = "character", help = "TSV with columns chrom, pos"),
  make_option("--out",       type = "character"),
  make_option("--window",    type = "double", default = 1e6,
              help = "Search radius for the nearest gene when none overlaps")
)))

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "gene_utils.R"))

pos <- read.delim(opt$positions, stringsAsFactors = FALSE, colClasses = c(chrom = "character"))
res <- do.call(rbind, lapply(seq_len(nrow(pos)), function(i) {
  g <- tryCatch(nearest_gene(pos$chrom[i], as.numeric(pos$pos[i]), opt$window),
                error = function(e) list(gene = ".", biotype = ".", distance = NA_integer_))
  data.frame(chrom = pos$chrom[i], pos = pos$pos[i],
             Gene = g$gene, Gene_Biotype = g$biotype, Gene_Distance_bp = g$distance,
             stringsAsFactors = FALSE)
}))

write.table(res, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("[gene_annotate] %d position(s); %d assigned to a named gene\n",
            nrow(res), sum(res$Gene != ".")))
