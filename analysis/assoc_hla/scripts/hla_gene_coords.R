#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose : GRCh38 coordinates for the HLA genes, from Ensembl rather than typed
#           in by hand. The marker positions in the .bim are derived from these,
#           so they have to be auditable and reproducible, not folklore.
#
#           TWO GENES ARE NOT HERE AND CANNOT BE. HLA-DRB3 and HLA-DRB4 exist
#           only on ALT haplotypes in GRCh38 -- there is no primary-assembly
#           locus to give them. Any position assigned to them would be invented.
#           They are excluded from the primary analysis for that reason among
#           others; see docs/METHODS.md.
# Component: assoc_hla
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse); library(EnsDb.Hsapiens.v86); library(AnnotationFilter)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option('--genes', type = 'character',
              help = 'comma-separated bare gene symbols, e.g. A,B,C,DRB1'),
  make_option('--out', type = 'character', help = 'output TSV'))))

bare <- strsplit(opt$genes, ',')[[1]]
want <- paste0('HLA-', bare)

g <- genes(EnsDb.Hsapiens.v86, filter = GeneNameFilter(want))
d <- data.frame(gene       = sub('^HLA-', '', mcols(g)$gene_name),
                chrom      = as.character(seqnames(g)),
                gene_start = start(g),
                gene_end   = end(g),
                stringsAsFactors = FALSE)

# Only the primary assembly. Ensembl also returns ALT-scaffold copies for some
# of these, and a marker placed on an ALT contig would not plot on chr6 and
# would not be excluded by SAIGE's --LOCO chr6 hold-out.
d <- d[d$chrom == '6', ]
d <- d[!duplicated(d$gene), ]
d <- d[order(d$gene_start), ]

missing <- setdiff(bare, d$gene)
if (length(missing))
  cat(sprintf('[hla_gene_coords] NOT on the GRCh38 primary assembly, dropped: %s\n',
              paste(missing, collapse = ', ')), file = stderr())

write.table(d, opt$out, sep = '\t', quote = FALSE, row.names = FALSE)
cat(sprintf('[hla_gene_coords] %d gene(s) -> %s (chr6 %s-%s)\n',
            nrow(d), opt$out, format(min(d$gene_start), big.mark = ','),
            format(max(d$gene_end), big.mark = ',')))
