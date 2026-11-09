#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose : Read back what a reused SAIGE null model actually contains, so the
#           random model is never run against a null fitted on other samples,
#           another phenotype or other covariates.
#
#           The null models are the genome-wide scan's (analysis/assoc_saige,
#           model 'main'); this component does not refit them. SAIGE step 2
#           matches genotypes to the null by sample ID and silently drops
#           samples it cannot match, so the sample set, the case count and the
#           covariate list are written out here and checked against this
#           component's phenotype file in VALIDATE_INPUTS.
#
# Output  : <out-prefix>.samples.tsv  sample_id, y (0/1) in the model's order
#           <out-prefix>.json         trait type, LOCO, n, cases, covariates,
#                                     variance ratio(s)
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) stop("missing ", flag)
    args[[i + 1]]
}
rda <- get_arg("--rda")
vr_file <- get_arg("--variance-ratio")
out <- get_arg("--out-prefix")

e <- new.env()
load(rda, envir = e)
if (!exists("modglmm", envir = e)) stop("no 'modglmm' object in ", rda)
m <- get("modglmm", envir = e)

ids <- as.character(m$sampleID)
y <- as.numeric(m$y)
if (length(ids) != length(y)) stop("sampleID and y differ in length in ", rda)
write.table(data.frame(sample_id = ids, y = y), paste0(out, ".samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# SAIGE stores the fixed-effect names after its QR step; the first is the
# intercept (named 'minus1' by SAIGE 1.1.9), the rest are the covariates.
coef_names <- names(m$coefficients)
# One row per variance-ratio category: value, type ('null' = the single
# non-categorical ratio), category index.
vr <- as.numeric(read.table(vr_file, header = FALSE, stringsAsFactors = FALSE)[[1]])

q <- function(s) paste0('"', s, '"')
json <- paste0(
    "{\n",
    '  "rda": ', q(normalizePath(rda)), ",\n",
    '  "trait_type": ', q(m$traitType), ",\n",
    '  "loco": ', tolower(as.character(isTRUE(m$LOCO))), ",\n",
    '  "n": ', length(ids), ",\n",
    '  "cases": ', sum(y == 1), ",\n",
    '  "controls": ', sum(y == 0), ",\n",
    '  "coefficients": [', paste(q(coef_names), collapse = ", "), "],\n",
    '  "variance_ratio": [', paste(format(vr, digits = 8), collapse = ", "), "]\n",
    "}\n")
writeLines(json, paste0(out, ".json"), sep = "")
cat(sprintf("[null_model_info] %s: n=%d cases=%d LOCO=%s covariates=%s\n",
            basename(rda), length(ids), sum(y == 1), isTRUE(m$LOCO),
            paste(coef_names[-1], collapse = ",")))
