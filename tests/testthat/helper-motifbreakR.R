## Shared fixtures for the motifbreakR tests. Variants are built from small
## BED files or from the bundled example results so that no SNPlocs lookup or
## network access is needed.

hg19 <- function() {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
}

## write BED lines (chrom, 0-based start, end, name) and import them
variants_from_bed <- function(lines, genome = hg19()) {
  bed <- withr::local_tempfile(fileext = ".bed", .local_envir = parent.frame())
  writeLines(paste(lines, 0, "+", sep = "\t"), bed)
  snps.from.file(bed, search.genome = genome, format = "bed")
}

## the variants of the bundled example.results, as motifbreakR input
example_variants <- function() {
  hg19()
  data("example.results", package = "motifbreakR", envir = environment())
  v <- example.results[!duplicated(paste(example.results$SNP_id, example.results$ALT))]
  mcols(v) <- mcols(v)[, c("SNP_id", "REF", "ALT")]
  strand(v) <- "*"
  names(v) <- paste(v$SNP_id, v$ALT, sep = ":")
  attributes(v)$genome.package <- "BSgenome.Hsapiens.UCSC.hg19"
  v
}

hocomoco_core <- function() {
  subset(MotifDb::MotifDb,
         dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C"))
}

## the vignette's analysis of the example variants, run once and shared
.cache <- new.env()
example_run <- function() {
  if (is.null(.cache$example_run)) {
    .cache$example_run <- run_mb(example_variants(), hocomoco_core(), filterp = TRUE,
                                 threshold = 1e-4, bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
  }
  .cache$example_run
}

run_mb <- function(snps, motifs, ..., BPPARAM = BiocParallel::SerialParam()) {
  motifbreakR(snps, motifs, method = "ic", BPPARAM = BPPARAM, ...)
}

## best score of a scoring matrix over every window (either strand) that
## contains position `pos` of `x`; a deliberately simple reference
## implementation for checking motifbreakR's vectorised scan
brute_force_max <- function(score_matrix, x, pos) {
  w <- ncol(score_matrix)
  rc <- c(A = "T", C = "G", G = "C", T = "A")
  starts <- max(1L, pos - w + 1L):min(pos, length(x) - w + 1L)
  max(vapply(starts, function(st) {
    win <- x[st:(st + w - 1L)]
    fwd <- sum(score_matrix[cbind(match(win, rownames(score_matrix)), seq_len(w))])
    rev <- sum(score_matrix[cbind(match(rev(rc[win]), rownames(score_matrix)), seq_len(w))])
    max(fwd, rev)
  }, numeric(1)))
}

## key identifying a result row
result_key <- function(x) paste(x$SNP_id, x$ALT, x$providerName)
