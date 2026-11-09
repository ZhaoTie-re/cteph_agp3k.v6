# The facts of the HIBAG DQB1 model the deck quotes (training samples, classifiers, SNPs, alleles, out-of-bag
# accuracy), written once from the pinned model file to assets/hibag/DQB1_model.json (the model never changes).
# usage: Rscript hibag_model_summary.R <model .RData> <HIBAG library> <out .json>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[2], .libPaths()))
suppressPackageStartupMessages(library(HIBAG))
o <- get(load(a[1]))[["DQB1"]]
m <- hlaModelFromObj(o)
s <- summary(m, show = FALSE)
acc <- s$info["accuracy", "Mean"]
nsnp <- s$info["num.snp", "Mean"]
one <- function(x) if (length(x)) as.character(x[1]) else ""
cat(sprintf('{"gene": "DQB1", "training_samples": %d, "classifiers": %d, "snps": %d, "snps_used": %d, "alleles": %d, "allele_names": [%s], "oob_accuracy_pct": %.4f, "snps_per_classifier": %.2f, "hibag_version": "%s", "platform": "%s", "information": "%s"}\n',
            as.integer(o$n.samp), as.integer(s$num.classifier), length(o$snp.id), as.integer(s$num.snp), length(m$hla.allele),
            paste0('"', m$hla.allele, '"', collapse = ", "), acc, nsnp, as.character(packageVersion("HIBAG")),
            one(m$appendix$platform), one(m$appendix$information)), file = a[3])
