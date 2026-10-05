## prepareVariants() builds the sequence contexts that every motif is scanned
## over; scoring requires all contexts in a run to share one width.

seq_chr <- function(x) unname(as.character(x))

check_contexts <- function(x, k = 20L) {
  genome <- hg19()
  p <- prepareVariants(x, genome, k)

  expect_length(unique(width(p$contextRef)), 1)
  expect_length(unique(width(p$contextAlt)), 1)
  expect_identical(unique(width(p$contextRef)), unique(width(p$contextAlt)))

  ## each allele sits exactly at its recorded coordinate
  expect_identical(seq_chr(subseq(p$contextRef, start(p$contextCoordRef), end(p$contextCoordRef))),
                   seq_chr(x$REF))
  expect_identical(seq_chr(subseq(p$contextAlt, start(p$contextCoordAlt), end(p$contextCoordAlt))),
                   seq_chr(x$ALT))

  ## the flanks are the genome sequence, at least k bases on each side
  left <- start(p$contextCoordAlt) - 1L
  right <- width(p$contextAlt) - end(p$contextCoordAlt)
  expect_true(all(c(left, right) >= k))
  up <- getSeq(genome, GRanges(seqnames(x), IRanges(end = start(x) - 1L, width = left)))
  down <- getSeq(genome, GRanges(seqnames(x), IRanges(start = end(x) + 1L, width = right)))
  expect_identical(seq_chr(subseq(p$contextAlt, 1L, left)), seq_chr(up))
  expect_identical(seq_chr(subseq(p$contextAlt, end(p$contextCoordAlt) + 1L, width(p$contextAlt))),
                   seq_chr(down))
}

test_that("contexts are consistent when the longest allele has even length", {
  ## regression: an even maximum allele length gave unequal context widths
  check_contexts(variants_from_bed(c("chr19\t38778914\t38778915\tchr19:38778915:A:AT",
                                     "chr19\t38778914\t38778916\tchr19:38778915:AC:A",
                                     "chr19\t38778914\t38778915\tchr19:38778915:A:G",
                                     "chr19\t38778914\t38778916\tchr19:38778915:AC:GT")))
})

test_that("contexts are consistent when the longest allele has odd length", {
  check_contexts(variants_from_bed(c("chr19\t38778914\t38778915\tchr19:38778915:A:ATG",
                                     "chr19\t38778914\t38778915\tchr19:38778915:A:G")))
})

test_that("contexts are consistent for SNVs only", {
  check_contexts(example_variants())
})
