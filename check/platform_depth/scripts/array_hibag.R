# DQB1 imputed from the SNP arrays with HIBAG and its published Asian two-field model (called by array_hibag.py).
# The model (assets/hibag/Asian-HLA4-hg38.RData: GSK HLARES of East and South Asian ancestry and HapMap CHB+JPT, HLA
# typed by classical methods; the hg38 file is the hg19 model lifted by position, so match.type = "Position" as its
# page advises) never saw our WGS typing: it is an independent measurement of DQB1 in cases and controls alike.
# usage: Rscript array_hibag.R <bfile prefix> <model .RData> <HIBAG library> <out prefix> <threads>
a <- commandArgs(trailingOnly = TRUE)
bfile <- a[1]; model_rdata <- a[2]; lib <- a[3]; out <- a[4]; threads <- as.integer(a[5])
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(HIBAG))
mobj <- get(load(model_rdata))
model <- hlaModelFromObj(mobj[["DQB1"]])
geno <- hlaBED2Geno(bed.fn = paste0(bfile, ".bed"), fam.fn = paste0(bfile, ".fam"), bim.fn = paste0(bfile, ".bim"),
                    assembly = "hg38")
n_found <- sum(model$snp.position %in% geno$snp.position)
pred <- hlaPredict(model, geno, type = "response+prob", match.type = "Position", cl = threads, verbose = FALSE)
write.table(pred$value, paste0(out, ".best.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pp <- pred$postprob
long <- do.call(rbind, lapply(seq_len(ncol(pp)), function(j) {
  k <- which(pp[, j] >= 0.001)
  data.frame(sample_id = colnames(pp)[j], pair = rownames(pp)[k], prob = pp[k, j])
}))
write.table(long, paste0(out, ".postprob.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(sprintf("hibag_version\t%s", as.character(packageVersion("HIBAG"))),
             sprintf("model_snps\t%d", length(model$snp.position)), sprintf("model_snps_on_array\t%d", n_found),
             sprintf("model_alleles\t%d", length(model$hla.allele)),
             sprintf("array_snps\t%d", length(geno$snp.position)), sprintf("samples\t%d", nrow(pred$value))),
           paste0(out, ".run.tsv"))
