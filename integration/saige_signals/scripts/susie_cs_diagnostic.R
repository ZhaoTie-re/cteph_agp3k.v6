#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose : Re-derive one peak's SuSiE fit exactly as analysis/_shared/scripts/
#           susie_finemap.R produced it, and list every 95 % credible set BEFORE
#           the purity filter. The pipeline keeps only sets with min |r| >= 0.5
#           and records nothing about the ones it dropped, so "no credible set"
#           cannot be explained from its outputs alone. Two mechanisms exist:
#             null fit  every single-effect prior variance V shrank to 0 -> all
#                       PIP are 0 and susie_get_cs returns nothing
#             impure    effects were fitted, but each 95 % set had min |r| < 0.5
#           This script separates them and gives the size and purity of the
#           dropped set(s), plus per-variant alpha so "outside the set" can be
#           stated as PIP below the set's cut. Several credible sets per peak
#           are supported throughout.
# Inputs  : the pipeline's own <peak>.sumstat.tsv, .ld_matrix.tsv, .ld_matrix.vars,
#           .pip.tsv, .cs.tsv, .susie.json; the fit is asserted to reproduce the
#           stored PIPs (default tolerance 1e-3; observed <= 1e-13).
#           --lead-id, --n, --L, --coverage and --min-abs-corr default to the
#           values recorded in --json, i.e. the arguments of the stored fit.
# Outputs : --out-effects  one row per effect l = 1..L
#           --out-variants one row per ID in --variants that lies in the fitted
#                          window (IDs outside it are counted, not fatal)
#           --out-fit      fit-level JSON
# Runtime : conda env r_work (R 4.4, susieR 0.14.2, optparse, jsonlite)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(susieR)
  library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstat",      type = "character"),
  make_option("--ld-matrix",    type = "character"),
  make_option("--ld-vars",      type = "character"),
  make_option("--pip",          type = "character", help = "the pipeline's <peak>.pip.tsv"),
  make_option("--cs",           type = "character", help = "the pipeline's <peak>.cs.tsv"),
  make_option("--json",         type = "character", help = "the pipeline's <peak>.susie.json"),
  make_option("--cohort",       type = "character"),
  make_option("--locus-id",     type = "character"),
  make_option("--lead-id",      type = "character", default = NULL, help = "default: lead_id in --json"),
  make_option("--n",            type = "double",    default = NULL, help = "default: n_gwas in --json"),
  make_option("--L",            type = "integer",   default = NULL, help = "default: L in --json"),
  make_option("--coverage",     type = "double",    default = NULL, help = "default: coverage in --json"),
  make_option("--min-abs-corr", type = "double",    default = NULL, help = "default: min_abs_corr in --json"),
  make_option("--variants",     type = "character", help = "one variant ID per line"),
  make_option("--tol",          type = "double",  default = 1e-3),
  make_option("--out-effects",  type = "character"),
  make_option("--out-variants", type = "character"),
  make_option("--out-fit",      type = "character")
)))

fail <- function(msg) { cat("[susie-diag] ", msg, "\n", sep = ""); quit(save = "no", status = 1) }
if (as.character(packageVersion("susieR")) != "0.14.2")
  fail(sprintf("susieR %s, the pipeline used 0.14.2", packageVersion("susieR")))

# the stored fit's own arguments are the defaults, so a refit cannot drift from it
meta <- jsonlite::fromJSON(opt$json)
if (is.null(opt$`lead-id`))      opt$`lead-id`      <- meta$lead_id
if (is.null(opt$n))              opt$n              <- as.numeric(meta$n_gwas)
if (is.null(opt$L))              opt$L              <- as.integer(meta$L)
if (is.null(opt$coverage))       opt$coverage       <- as.numeric(meta$coverage)
if (is.null(opt$`min-abs-corr`)) opt$`min-abs-corr` <- as.numeric(meta$min_abs_corr)

# ── the pipeline's preprocessing, verbatim ───────────────────────────────────
ss <- read.delim(opt$sumstat, stringsAsFactors = FALSE, check.names = FALSE)
names(ss)[names(ss) == "LOG(OR)_SE"] <- "SE"
if (!"SE" %in% names(ss)) fail("summary statistics carry no standard-error column")

ld  <- as.matrix(read.table(opt$`ld-matrix`, header = FALSE))
vid <- read.csv(opt$`ld-vars`, stringsAsFactors = FALSE)$ID
if (nrow(ld) != length(vid)) fail(sprintf("LD is %dx%d but %d variant IDs", nrow(ld), ncol(ld), length(vid)))
rownames(ld) <- colnames(ld) <- vid

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
R[!is.finite(R)] <- 0
diag(R) <- 1
R <- (R + t(R)) / 2

fit <- susie_rss(bhat = bhat, shat = shat, R = R, n = opt$n,
                 L = opt$L, estimate_residual_variance = FALSE)
if (!isTRUE(fit$converged)) fail("susie_rss did not converge")

# ── reproduce the stored PIP and the stored (post-purity) credible set ───────
pip_new <- as.numeric(susie_get_pip(fit))
pip_old <- read.delim(opt$pip, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(pip_old$ID, ss$ID)) fail("variant order differs from the stored pip.tsv")
max_dpip <- max(abs(pip_new - as.numeric(pip_old$pip)))
if (!is.finite(max_dpip) || max_dpip > opt$tol)
  fail(sprintf("stored PIPs not reproduced: max |delta| = %g", max_dpip))
if (!identical(as.logical(meta$converged), TRUE)) fail("stored fit is recorded as not converged")

set.seed(1)
cs_all  <- susie_get_cs(fit, Xcorr = R, coverage = opt$coverage, min_abs_corr = 0, n_purity = Inf)
cs_pure <- susie_get_cs(fit, Xcorr = R, coverage = opt$coverage, min_abs_corr = opt$`min-abs-corr`)
n_pure  <- length(cs_pure$cs)
if (n_pure != as.integer(meta$n_credible_sets))
  fail(sprintf("%d pure set(s) now, %d stored", n_pure, meta$n_credible_sets))

# every stored set (there may be several per peak) must be reproduced by one
# surviving set, with the same coverage and purity
stored_cs <- read.delim(opt$cs, stringsAsFactors = FALSE, check.names = FALSE)
stored_sets <- lapply(seq_len(nrow(stored_cs)),
                      function(i) strsplit(stored_cs$variants[i], ",", fixed = TRUE)[[1]])
if (nrow(stored_cs) != n_pure) fail(sprintf("%d stored set(s) but %d survive the refit", nrow(stored_cs), n_pure))
for (i in seq_along(stored_sets)) {
  hit <- which(vapply(cs_pure$cs, function(idx) setequal(ss$ID[idx], stored_sets[[i]]), logical(1)))
  if (length(hit) != 1) fail(sprintf("stored set %d is not reproduced by exactly one surviving set", i))
  if (abs(cs_pure$coverage[hit] - stored_cs$coverage[i]) > 1e-6) fail("stored coverage not reproduced")
  if (abs(cs_pure$purity$min.abs.corr[hit] - stored_cs$purity_min_abs_corr[i]) > 1e-6)
    fail("stored purity not reproduced")
}
in_stored <- function(idx) any(vapply(stored_sets, function(m) setequal(ss$ID[idx], m), logical(1)))

# ── every effect's own 95 % set (susieR's in_CS rule), exact purity ──────────
alpha <- fit$alpha                      # L x p
V     <- as.numeric(fit$V)
active <- V > 1e-9                      # susie_get_pip / susie_get_cs include_idx
L <- nrow(alpha)
own_set <- function(a) {
  o <- order(a, decreasing = TRUE)
  n <- sum(cumsum(a[o]) < opt$coverage) + 1
  o[seq_len(min(n, length(o)))]
}
purity_of <- function(idx) {
  if (length(idx) == 1) return(c(1, 1, 1))
  m <- abs(R[idx, idx, drop = FALSE])
  v <- m[upper.tri(m)]
  c(min(v), mean(v), median(v))
}
sets <- lapply(seq_len(L), function(l) if (active[l]) own_set(alpha[l, ]) else integer())
# effects with an identical set are one set after susie_get_cs()'s dedup
dedup_of <- rep(NA_integer_, L)
for (l in seq_len(L)) if (active[l]) for (k in seq_len(l - 1))
  if (active[k] && is.na(dedup_of[l]) && setequal(sets[[l]], sets[[k]])) dedup_of[l] <- k
# cross-check against susie_get_cs(min_abs_corr = 0): every returned set is one of ours
for (i in seq_along(cs_all$cs)) {
  hit <- which(vapply(sets, function(s) length(s) && setequal(s, cs_all$cs[[i]]), logical(1)))
  if (!length(hit)) fail(sprintf("set %d of susie_get_cs is not any effect's own set", i))
  p <- purity_of(cs_all$cs[[i]])
  if (abs(p[1] - cs_all$purity$min.abs.corr[i]) > 1e-6) fail("purity recomputation disagrees with susieR")
}
pure_effects <- if (n_pure) as.integer(cs_pure$cs_index) else integer()

table_ids <- readLines(opt$variants)
table_ids <- table_ids[nzchar(table_ids)]
# a table variant can lie outside this cohort's fitted window (it entered the
# table through another cohort's credible set); it is reported upstream as
# "not fine-mapped", so here it is only counted
n_outside <- sum(!table_ids %in% ss$ID)
if (n_outside) cat(sprintf("[susie-diag] %d table variant(s) outside the fitted window\n", n_outside))
table_ids <- table_ids[table_ids %in% ss$ID]
lead_idx  <- match(opt$`lead-id`, ss$ID)

rows <- lapply(seq_len(L), function(l) {
  idx <- sets[[l]]
  p   <- if (length(idx)) purity_of(idx) else c(NA, NA, NA)
  in_table <- if (length(idx)) intersect(ss$ID[idx], table_ids) else character()
  data.frame(cohort = opt$cohort, locus_id = opt$`locus-id`, lead_id = opt$`lead-id`,
             n_gwas = opt$n, n_variants = length(ss$ID), L = L, niter = fit$niter,
             converged = isTRUE(fit$converged),
             n_effects = sum(active), n_sets_pre = length(cs_all$cs), n_sets_pure = n_pure,
             n_credible_sets_stored = as.integer(meta$n_credible_sets), max_abs_dpip = max_dpip,
             effect = l, V = V[l], active = active[l], dedup_of = dedup_of[l],
             set_size = if (length(idx)) length(idx) else NA_integer_,
             set_coverage = if (length(idx)) sum(alpha[l, idx]) else NA_real_,
             min_abs_corr = p[1], mean_abs_corr = p[2], median_abs_corr = p[3],
             purity_pass = if (length(idx)) p[1] >= opt$`min-abs-corr` else NA,
             stored = if (length(idx)) in_stored(idx) else FALSE,
             contains_lead = if (length(idx)) lead_idx %in% idx else FALSE,
             top_variant = if (length(idx)) ss$ID[idx[which.max(alpha[l, idx])]] else NA_character_,
             top_alpha = if (length(idx)) max(alpha[l, idx]) else NA_real_,
             n_table_variants = length(in_table),
             table_variants = paste(in_table, collapse = ","),
             stringsAsFactors = FALSE)
})
eff <- do.call(rbind, rows)
write.table(eff, opt$`out-effects`, sep = "\t", quote = FALSE, row.names = FALSE)

# ── per-variant view for the table's variants ────────────────────────────────
pos <- match(table_ids, ss$ID)
# surviving sets, numbered as susie_get_cs returns them (1..n_pure); each has an
# effect (cs_index) and a cut = the smallest alpha among its members
set_cut <- if (n_pure) vapply(seq_len(n_pure), function(k) min(alpha[pure_effects[k], cs_pure$cs[[k]]]), numeric(1)) else numeric()
stored_index_of <- function(id) {              # which stored set (row of cs.tsv) holds the variant
  k <- which(vapply(stored_sets, function(m) id %in% m, logical(1)))
  if (length(k)) k[1] else NA_integer_
}
vrows <- lapply(seq_along(table_ids), function(i) {
  j <- pos[i]
  pre <- which(vapply(sets, function(s) j %in% s, logical(1)))
  member_of <- if (n_pure) which(vapply(cs_pure$cs, function(idx) j %in% idx, logical(1))) else integer()
  # the surviving set this variant is closest to: the one whose effect gives it
  # the largest alpha (its own set when it is a member)
  k <- if (length(member_of)) member_of[1] else if (n_pure) which.max(vapply(seq_len(n_pure), function(k) alpha[pure_effects[k], j], numeric(1))) else NA_integer_
  out <- data.frame(cohort = opt$cohort, locus_id = opt$`locus-id`, ID = table_ids[i],
                    is_lead = table_ids[i] == opt$`lead-id`,
                    pip_stored = as.numeric(pip_old$pip[j]), pip_new = pip_new[j],
                    stored_cs = stored_index_of(table_ids[i]),
                    stringsAsFactors = FALSE)
  for (l in seq_len(L)) out[[sprintf("alpha_%d", l)]] <- alpha[l, j]
  out$pre_sets <- paste(pre, collapse = ";")
  out$pure_sets <- paste(member_of, collapse = ";")
  out$pure_set <- if (length(member_of)) pure_effects[member_of[1]] else NA_integer_
  out$nearest_set <- if (is.na(k)) NA_integer_ else as.integer(k)
  out$alpha_pure <- if (is.na(k)) NA_real_ else alpha[pure_effects[k], j]
  out$alpha_cut_pure <- if (is.na(k)) NA_real_ else set_cut[k]
  out
})
vtab <- if (length(vrows)) do.call(rbind, vrows) else
  data.frame(cohort = character(), locus_id = character(), ID = character(), is_lead = logical(),
             pip_stored = numeric(), pip_new = numeric(), stored_cs = integer(),
             setNames(as.data.frame(matrix(numeric(), 0, L)), sprintf("alpha_%d", seq_len(L))),
             pre_sets = character(), pure_sets = character(), pure_set = integer(), nearest_set = integer(),
             alpha_pure = numeric(), alpha_cut_pure = numeric())
write.table(vtab, opt$`out-variants`, sep = "\t", quote = FALSE, row.names = FALSE)

fitmeta <- list(cohort = opt$cohort, locus_id = opt$`locus-id`, lead_id = opt$`lead-id`,
                n_variants = length(ss$ID), n_gwas = opt$n, L = L, coverage = opt$coverage,
                min_abs_corr = opt$`min-abs-corr`, n_table_outside = n_outside, niter = fit$niter,
                converged = isTRUE(fit$converged), V = V, sigma2 = fit$sigma2,
                elbo_last = tail(fit$elbo, 1), n_effects = sum(active),
                n_sets_pre = length(cs_all$cs), n_sets_pure = n_pure,
                n_credible_sets_stored = as.integer(meta$n_credible_sets),
                max_abs_dpip = max_dpip, susieR_version = as.character(packageVersion("susieR")),
                purity = "exact (n_purity = Inf, all pairs)")
writeLines(jsonlite::toJSON(fitmeta, auto_unbox = TRUE, pretty = TRUE, digits = NA), opt$`out-fit`)
cat(sprintf("[susie-diag] %s %s: %d variants, %d effect(s) fitted, %d set(s) before purity, %d kept, max |dPIP| = %.1e\n",
            opt$cohort, opt$`locus-id`, length(ss$ID), sum(active), length(cs_all$cs), n_pure, max_dpip))
