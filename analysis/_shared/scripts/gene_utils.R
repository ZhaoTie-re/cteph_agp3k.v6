#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose : Shared gene-model helpers for this component. Sourced by
#           region_tracks.R (the regional-plot gene track) and gene_annotate.R
#           (the `Gene` column of the lead-variant table), so both apply exactly
#           the same notion of "a gene worth naming".
# Runtime : conda env `r_work` — EnsDb.Hsapiens.v86 (GRCh38, Ensembl 86)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
  library(GenomicRanges)
})

EDB <- EnsDb.Hsapiens.v86

# Clone-based and unnamed models. Ensembl carries many gene entries whose
# "name" is only the sequencing clone or scaffold accession they were called
# from — AC007347.1, RP11-357N13.3, CTD-2201I18.1. They convey nothing to a
# reader and in a dense window can outnumber the named genes several to one. They are
# suppressed everywhere; HGNC-approved non-coding symbols (LINC01234, MIR548,
# SNORD116) are NOT clone names and are kept.
CLONE_PATTERNS <- c(
  "^(AC|AL|AP|AF|BX|CR|CU|FP|FO|LL|Z)[0-9]{6}\\.[0-9]+$",  # clone accessions
  "^(RP|CTD|CTA|CTB|CTC|LA|XX|KB|DK|WI)[0-9]*-",           # library-plate clones
  "^Clone_", "^U[0-9]+\\.[0-9]+$"
)


# Biotypes that carry no interpretable signal on a regional plot even when the
# model happens to have a symbol.
DROP_BIOTYPES <- c("TEC", "artifact")

#' Genes in a window that are worth naming.
#'
#' The rule is: the gene must carry an OFFICIAL symbol — not a clone accession,
#' not a bare ENSG id — and must not be of a nuisance biotype. Protein-coding
#' genes and HGNC-approved non-coding genes (LINC*, MIR*, SNOR*) both qualify;
#' AC007347.1 and RP11-357N13.3 do not.
#'
#' @param chrom  chromosome, with or without a 'chr' prefix
#' @param start,end  window in bp
#' @param drop_pseudo  also drop *_pseudogene biotypes (default TRUE)
informative_genes <- function(chrom, start, end, drop_pseudo = TRUE) {
  chr <- sub("^chr", "", as.character(chrom))
  gr <- ensembldb::genes(
    EDB,
    filter = AnnotationFilterList(
      SeqNameFilter(chr),
      GeneStartFilter(end, condition = "<="),
      GeneEndFilter(start, condition = ">=")))
  if (length(gr) == 0) return(gr)

  nm <- gr$gene_name
  keep <- !is.na(nm) & nzchar(nm) & !grepl("^ENSG[0-9]+", nm)
  for (p in CLONE_PATTERNS) keep <- keep & !grepl(p, nm)
  keep <- keep & !(gr$gene_biotype %in% DROP_BIOTYPES)
  if (drop_pseudo) keep <- keep & !grepl("pseudogene$", gr$gene_biotype)
  gr[keep]
}


#' Exon structure of one representative transcript per gene.
#'
#' The representative is the transcript with the greatest total exonic length,
#' preferring `protein_coding` when the gene has one — a gene drawn from its
#' longest coding transcript shows the structure a reader expects, whereas
#' drawing every transcript of a multi-transcript gene would need one row each.
#'
#' @return data.frame: gene, tx_id, tx_biotype, strand, gene_start, gene_end,
#'   exon_start, exon_end, exon_rank  (one row per exon)
representative_exons <- function(gr_genes) {
  if (length(gr_genes) == 0) {
    return(data.frame(gene = character(), tx_id = character(), tx_biotype = character(),
                      strand = character(), gene_start = numeric(), gene_end = numeric(),
                      exon_start = numeric(), exon_end = numeric(), exon_rank = integer(),
                      stringsAsFactors = FALSE))
  }
  ids <- gr_genes$gene_id
  tx <- ensembldb::transcripts(EDB, filter = GeneIdFilter(ids))
  if (length(tx) == 0) return(representative_exons(gr_genes[0]))
  ex <- ensembldb::exonsBy(EDB, by = "tx", filter = GeneIdFilter(ids))

  rows <- list()
  for (i in seq_along(gr_genes)) {
    gid  <- gr_genes$gene_id[i]
    gnm  <- gr_genes$gene_name[i]
    cand <- tx[tx$gene_id == gid]
    if (length(cand) == 0) next
    coding <- cand[cand$tx_biotype == "protein_coding"]
    if (length(coding) > 0) cand <- coding
    exlen <- vapply(cand$tx_id, function(t) {
      if (is.null(ex[[t]])) return(0) else sum(width(ex[[t]]))
    }, numeric(1))
    best <- cand$tx_id[which.max(exlen)]
    e <- ex[[best]]
    if (is.null(e) || length(e) == 0) next
    rows[[length(rows) + 1]] <- data.frame(
      gene = gnm, tx_id = best,
      tx_biotype = cand$tx_biotype[cand$tx_id == best][1],
      strand = as.character(strand(gr_genes)[i]),
      gene_start = start(gr_genes)[i], gene_end = end(gr_genes)[i],
      exon_start = start(e), exon_end = end(e), exon_rank = seq_along(e),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(representative_exons(gr_genes[0]))
  do.call(rbind, rows)
}


#' Gene assignment for a single position.
#'
#' Overlapping gene if there is one (preferring protein_coding, then the longest
#' model), otherwise the nearest. `distance` is 0 when the variant is inside the
#' gene; otherwise it is the gap in bp, NEGATIVE when the variant lies before the
#' gene's start and POSITIVE when it lies after the gene's end.
nearest_gene <- function(chrom, pos, window = 1e6) {
  gr <- informative_genes(chrom, max(1, pos - window), pos + window)
  if (length(gr) == 0) return(list(gene = ".", biotype = ".", distance = NA_integer_))
  s <- start(gr); e <- end(gr)
  inside <- which(pos >= s & pos <= e)
  if (length(inside)) {
    pc <- inside[gr$gene_biotype[inside] == "protein_coding"]
    k <- if (length(pc)) pc[which.max(e[pc] - s[pc])] else inside[which.max(e[inside] - s[inside])]
    return(list(gene = gr$gene_name[k], biotype = gr$gene_biotype[k], distance = 0L))
  }
  d <- ifelse(pos < s, pos - s, pos - e)
  k <- which.min(abs(d))
  list(gene = gr$gene_name[k], biotype = gr$gene_biotype[k], distance = as.integer(d[k]))
}
