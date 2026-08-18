#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose : SuSiE fine-mapping of one genome-wide locus from summary statistics
#           plus an in-sample LD matrix. Emits per-variant PIP and the 95 %
#           credible sets with their purity.
# Runtime : conda env `r_work` (R 4.4.3, susieR 0.14.2, optparse 1.8.2)
# Used by : assoc_plink2.nf  process SUSIE
#
# Four corrections over the reference implementation, each a real defect:
#   1. susie_get_cs(..., Xcorr = R), NOT X = R. In susieR's signature `X` is the
#      genotype matrix and `Xcorr` the correlation matrix; passing an LD matrix
#      as X makes susieR read p variants as p "samples", so the purity it
#      reports is not credible-set purity and min_abs_corr never filters.
#   2. Align bhat/shat to the LD matrix BY VARIANT ID, with an assertion. The
#      reference relied on the extract list and the .bim happening to agree in
#      order and never checked; a single mismatch silently permutes R.
#   3. n = the GWAS analysis sample size, not the LD reference size. The
#      reference passed the number of LD samples, which distorts the implied
#      z-to-effect scaling.
#   4. LD from all analysed samples, not controls only: R must describe the
#      sample the summary statistics came from.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(susieR)
  library(Matrix)
  library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstat",   type = "character"),
  make_option("--ld-matrix", type = "character"),
  make_option("--ld-vars",   type = "character"),
  make_option("--cohort",    type = "character"),
  make_option("--locus-id",  type = "character"),
  make_option("--lead-id",   type = "character"),
  make_option("--n",         type = "double", help = "GWAS analysis N (see correction 3)"),
  make_option("--L",         type = "integer", default = 10L),
  make_option("--coverage",  type = "double",  default = 0.95),
  make_option("--min-abs-corr", type = "double", default = 0.5),
  make_option("--out-pip",   type = "character"),
  make_option("--out-cs",    type = "character"),
  make_option("--out-json",  type = "character")
)))

fail <- function(msg) { cat("[susie] ", msg, "\n", sep = ""); quit(save = "no", status = 1) }

ss <- read.delim(opt$sumstat, stringsAsFactors = FALSE, check.names = FALSE)
names(ss)[names(ss) == "LOG(OR)_SE"] <- "SE"
if (!"SE" %in% names(ss)) fail("summary statistics carry no standard-error column")

ld  <- as.matrix(read.table(opt$`ld-matrix`, header = FALSE))
vid <- read.csv(opt$`ld-vars`, stringsAsFactors = FALSE)$ID
if (nrow(ld) != length(vid)) fail(sprintf("LD is %dx%d but %d variant IDs", nrow(ld), ncol(ld), length(vid)))
rownames(ld) <- colnames(ld) <- vid

# ── correction 2: align by ID, then assert ─────────────────────────────────
common <- intersect(ss$ID, vid)
if (length(common) < 2) fail(sprintf("only %d variants shared between sumstats and LD", length(common)))
ss <- ss[match(common, ss$ID), , drop = FALSE]
R  <- ld[common, common, drop = FALSE]
stopifnot(identical(ss$ID, rownames(R)), identical(ss$ID, colnames(R)))

bhat <- log(as.numeric(ss$OR))
shat <- as.numeric(ss$SE)
keep <- is.finite(bhat) & is.finite(shat) & shat > 0
if (sum(keep) < 2) fail("fewer than 2 variants with usable effect estimates")
ss <- ss[keep, , drop = FALSE]; bhat <- bhat[keep]; shat <- shat[keep]
R <- R[keep, keep, drop = FALSE]

# susie_rss requires a symmetric, finite R. PLINK writes 'nan' where a variant
# is monomorphic in the LD sample; those rows carry no information, so they are
# zeroed off-diagonal rather than allowed to poison the decomposition.
R[!is.finite(R)] <- 0
diag(R) <- 1
R <- (R + t(R)) / 2

fit <- tryCatch(
  susie_rss(bhat = bhat, shat = shat, R = R, n = opt$n,
            L = opt$L, estimate_residual_variance = FALSE),
  error = function(e) { cat("[susie] susie_rss failed: ", conditionMessage(e), "\n", sep = ""); NULL })

pip_tbl <- data.frame(cohort = opt$cohort, locus_id = opt$`locus-id`,
                      ID = ss$ID, CHROM = ss$CHROM, POS = ss$POS,
                      A1 = if ("A1" %in% names(ss)) ss$A1 else NA,
                      A1_FREQ = if ("A1_FREQ" %in% names(ss)) ss$A1_FREQ else NA,
                      P = ss$P, OR = ss$OR, SE = shat,
                      is_lead = ss$ID == opt$`lead-id`,
                      pip = NA_real_, cs = NA_integer_,
                      stringsAsFactors = FALSE)
cs_tbl <- data.frame()

if (!is.null(fit)) {
  pip_tbl$pip <- as.numeric(susie_get_pip(fit))
  # ── correction 1: Xcorr =, not X = ───────────────────────────────────────
  cs <- susie_get_cs(fit, Xcorr = R, coverage = opt$coverage,
                     min_abs_corr = opt$`min-abs-corr`)
  if (length(cs$cs) > 0) {
    rows <- lapply(seq_along(cs$cs), function(i) {
      idx <- cs$cs[[i]]
      pip_tbl$cs[idx] <<- i
      ordr <- idx[order(pip_tbl$pip[idx], decreasing = TRUE)]
      data.frame(cohort = opt$cohort, locus_id = opt$`locus-id`, cs = i,
                 size = length(idx),
                 coverage = cs$coverage[i],
                 purity_min_abs_corr = cs$purity$min.abs.corr[i],
                 purity_mean_abs_corr = cs$purity$mean.abs.corr[i],
                 top_variant = pip_tbl$ID[ordr[1]],
                 top_pip = pip_tbl$pip[ordr[1]],
                 contains_lead = opt$`lead-id` %in% pip_tbl$ID[idx],
                 variants = paste(pip_tbl$ID[ordr], collapse = ","),
                 stringsAsFactors = FALSE)
    })
    cs_tbl <- do.call(rbind, rows)
  }
}

write.table(pip_tbl, opt$`out-pip`, sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(cs_tbl) == 0) {
  cs_tbl <- data.frame(cohort = character(), locus_id = character(), cs = integer(),
                       size = integer(), coverage = numeric(),
                       purity_min_abs_corr = numeric(), purity_mean_abs_corr = numeric(),
                       top_variant = character(), top_pip = numeric(),
                       contains_lead = logical(), variants = character())
}
write.table(cs_tbl, opt$`out-cs`, sep = "\t", quote = FALSE, row.names = FALSE)

meta <- list(cohort = opt$cohort, locus_id = opt$`locus-id`, lead_id = opt$`lead-id`,
             n_gwas = opt$n, n_variants = nrow(pip_tbl), L = opt$L,
             coverage = opt$coverage, min_abs_corr = opt$`min-abs-corr`,
             converged = if (is.null(fit)) FALSE else isTRUE(fit$converged),
             n_credible_sets = nrow(cs_tbl),
             susieR_version = as.character(packageVersion("susieR")))
writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), opt$`out-json`)

cat(sprintf("[susie] %s %s: %d variants, %d credible set(s), converged=%s\n",
            opt$cohort, opt$`locus-id`, nrow(pip_tbl), nrow(cs_tbl), meta$converged))
