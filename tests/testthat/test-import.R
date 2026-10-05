test_that("snps.from.file reads SNVs, indels and MNVs from BED", {
  x <- variants_from_bed(c("chr19\t38778914\t38778915\tchr19:38778915:A:AT",
                           "chr19\t38778914\t38778916\tchr19:38778915:AC:A",
                           "chr19\t38778914\t38778915\tchr19:38778915:A:G",
                           "chr19\t38778914\t38778916\tchr19:38778915:AC:GT"))

  expect_s4_class(x, "GRanges")
  expect_length(x, 4)
  expect_named(mcols(x), c("SNP_id", "REF", "ALT"))
  expect_s4_class(x$REF, "DNAStringSet")
  expect_setequal(paste(x$REF, x$ALT, sep = ">"), c("A>AT", "AC>A", "A>G", "AC>GT"))
  expect_identical(attributes(x)$genome.package, "BSgenome.Hsapiens.UCSC.hg19")
})

test_that("snps.from.file splits comma separated alternate alleles", {
  x <- variants_from_bed("chr19\t38778914\t38778915\tchr19:38778915:A:C,G,T")
  expect_length(x, 3)
  expect_setequal(as.character(x$ALT), c("C", "G", "T"))
})

test_that("snps.from.file reads indels from VCF and drops symbolic alleles", {
  vcf <- system.file("extdata", "chek2.vcf.gz", package = "motifbreakR")
  expect_warning(x <- snps.from.file(vcf, search.genome = hg19(), format = "vcf"),
                 "non-standard nucleotide codes")

  expect_gt(sum(width(x$REF) > 1 | width(x$ALT) > 1), 0)
  expect_false(any(grepl("[<>]", c(as.character(x$REF), as.character(x$ALT)))))
  expect_named(mcols(x), c("SNP_id", "REF", "ALT"))
})

test_that("variants.from.file is identical to snps.from.file", {
  vcf <- system.file("extdata", "chek2.vcf.gz", package = "motifbreakR")
  a <- suppressWarnings(snps.from.file(vcf, search.genome = hg19(), format = "vcf"))
  b <- suppressWarnings(variants.from.file(vcf, search.genome = hg19(), format = "vcf"))
  expect_identical(a, b)
})
