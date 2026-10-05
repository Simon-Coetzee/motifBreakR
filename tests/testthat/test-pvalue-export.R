test_that("calculatePvalue fills in p-values using the stored background", {
  data("example.results", package = "motifbreakR", envir = environment())
  x <- example.results[1:2]
  res <- calculatePvalue(x, granularity = 1e-3, BPPARAM = BiocParallel::SerialParam())

  expect_false(anyNA(res$pValueRef))
  expect_false(anyNA(res$pValueAlt))
  expect_true(all(res$pValueRef > 0 & res$pValueRef <= 1))
  expect_true(all(res$pValueEffect %in% c("strong", "weak")))
  expect_identical(res$scoreRef, x$scoreRef)
})

test_that("calculatePvalue requires results run with filterp = TRUE", {
  data("example.results", package = "motifbreakR", envir = environment())
  x <- example.results[1:2]
  mcols(x)$pValueRef <- NULL
  expect_error(calculatePvalue(x), "filterp=TRUE")
})

test_that("exportMBbed best_sig uses the lower p-value", {
  ## regression: best_sig used the p-value of the weaker allele
  load(system.file("extdata", "example.pvalue.rda", package = "motifbreakR"))
  bed <- withr::local_tempfile(fileext = ".bed")
  exportMBbed(rs1006140, file = bed, color = "best_sig")
  out <- rtracklayer::import(bed)
  expect_length(out, length(rs1006140))
  expect_equal(out$score, -log10(pmin(rs1006140$pValueRef, rs1006140$pValueAlt)),
               tolerance = 1e-6)
})

test_that("exportMBbed needs p-values for the significance colours", {
  data("example.results", package = "motifbreakR", envir = environment())
  bed <- withr::local_tempfile(fileext = ".bed")
  expect_error(exportMBbed(example.results, file = bed, color = "ref_sig"), "calculatePvalue")
  expect_no_error(exportMBbed(example.results, file = bed, color = "effect_size"))
})

test_that("exportMBtable writes one row per result", {
  data("example.results", package = "motifbreakR", envir = environment())
  for (fmt in c("tsv", "csv")) {
    f <- withr::local_tempfile(fileext = paste0(".", fmt))
    exportMBtable(example.results, file = f, format = fmt)
    tab <- utils::read.table(f, header = TRUE, sep = if (fmt == "tsv") "\t" else ",",
                             quote = "\"", comment.char = "")
    expect_equal(nrow(tab), length(example.results))
  }
})

test_that("plotMB draws without error", {
  skip_on_cran()
  skip_if_offline()
  hg19()
  data("example.results", package = "motifbreakR", envir = environment())
  f <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plotMB(example.results, rsid = "rs1006140", effect = "strong", altAllele = "C"))
})
