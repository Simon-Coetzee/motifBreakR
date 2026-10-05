result_columns <- c("SNP_id", "REF", "ALT", "varType",
                    "windowIdxRef", "motifStrandRef", "windowIdxAlt", "motifStrandAlt",
                    "motifID", "geneSymbol", "dataSource", "providerName", "providerId",
                    "seqMatch", "pwmConsensus", "pctRef", "pctAlt", "scoreRef", "scoreAlt",
                    "strongerIn", "alleleDiff", "alleleEffectSize", "effect")

test_that("motifbreakR reproduces the bundled example results", {
  data("example.results", package = "motifbreakR", envir = environment())
  res <- example_run()

  expect_setequal(result_key(res), result_key(example.results))
  m <- match(result_key(res), result_key(example.results))
  expect_equal(res$scoreRef, example.results$scoreRef[m])
  expect_equal(res$scoreAlt, example.results$scoreAlt[m])
  expect_identical(as.character(res$effect), as.character(example.results$effect[m]))
})

test_that("allele scores are the best match over windows containing the variant", {
  ## regression: windows not overlapping the variant used to be scored
  genome <- hg19()
  motifs <- hocomoco_core()
  res <- example_run()
  score_mats <- attributes(res)$scoremotifs
  ids <- paste(mcols(motifs)$providerId, mcols(motifs)$providerName)

  for (i in seq_along(res)) {
    r <- res[i]
    mat <- score_mats[[names(motifs)[match(paste(r$providerId, r$providerName), ids)]]]
    w <- ncol(mat)
    ref <- strsplit(as.character(getSeq(genome, as.character(seqnames(r)),
                                        start(r) - w, start(r) + w)), "")[[1]]
    alt <- ref
    alt[w + 1L] <- as.character(r$ALT)
    expect_equal(unname(r$scoreRef), brute_force_max(mat, ref, w + 1L), info = result_key(r))
    expect_equal(unname(r$scoreAlt), brute_force_max(mat, alt, w + 1L), info = result_key(r))
  }
})

test_that("filterp controls the p-value columns", {
  snps <- example_variants()
  motifs <- hocomoco_core()

  with_p <- example_run()
  expect_true(all(result_columns %in% names(mcols(with_p))))
  expect_true(all(c("pValueRef", "pValueAlt") %in% names(mcols(with_p))))
  expect_true(all(is.na(with_p$pValueRef)))

  ## regression: filterp = FALSE (the default) used to error
  without_p <- run_mb(snps, motifs, filterp = FALSE, threshold = 0.85)
  expect_true(all(result_columns %in% names(mcols(without_p))))
  expect_false(any(c("pValueRef", "pValueAlt") %in% names(mcols(without_p))))
  expect_true(all(pmax(without_p$pctRef, without_p$pctAlt) > 0.85))
})

test_that("results are well formed", {
  res <- example_run()
  expect_true(all(res$strongerIn == ifelse(res$alleleDiff > 0, "Alt", "Ref")))
  expect_true(all(res$motifStrandRef %in% c("+", "-")))
  expect_true(all(as.character(res$effect) %in% c("strong", "weak")))
  expect_true(all(width(res$seqMatch) == width(res$pwmConsensus)))
  expect_true(all(c("genome.package", "motifs", "scoremotifs", "bkg") %in% names(attributes(res))))
})

test_that("motifs sharing a providerId are scored separately", {
  ## regression: rows were selected by providerId, which is not unique
  data("hocomoco", package = "motifbreakR", envir = environment())
  motifs <- hocomoco[8:11]
  dup <- mcols(motifs)$providerId[2]
  expect_identical(sum(mcols(motifs)$providerId == dup), 2L)

  ## threshold 0 with show.neutral reports every motif for every variant
  snps <- example_variants()
  res <- run_mb(snps, motifs, threshold = 0, show.neutral = TRUE)
  expect_length(res, length(snps) * length(motifs))
  expect_false(anyDuplicated(result_key(res)) > 0)
  expect_setequal(res$providerName[res$providerId == dup],
                  mcols(motifs)$providerName[mcols(motifs)$providerId == dup])
})

test_that("indels and MNVs are scored", {
  snps <- variants_from_bed(c("chr19\t38778914\t38778915\tchr19:38778915:A:AT",
                              "chr19\t38778914\t38778916\tchr19:38778915:AC:A",
                              "chr19\t38778914\t38778915\tchr19:38778915:A:G",
                              "chr19\t38778914\t38778916\tchr19:38778915:AC:GT"))
  res <- run_mb(snps, subset(MotifDb::MotifDb, dataSource == "HOCOMOCOv11-core-A"), threshold = 0.8)
  expect_setequal(unique(res$varType), c("Insertion", "Deletion", "SNV", "Other"))
  expect_setequal(unique(res$SNP_id), snps$SNP_id)
})

test_that("a run with no hits returns NULL with a warning", {
  snps <- example_variants()[1]
  expect_warning(res <- run_mb(snps, hocomoco_core()[1:2], threshold = 0.999),
                 "No SNP/Motif Interactions reached threshold")
  expect_null(res)
})

test_that("bkg accepts frequencies and the named presets", {
  snps <- example_variants()
  motifs <- hocomoco_core()[1:40]
  for (b in list(c(A = 0.3, C = 0.2, G = 0.2, T = 0.3), "aggregate", "pwm")) {
    res <- run_mb(snps, motifs, threshold = 0.8, bkg = b)
    expect_equal(sum(attributes(res)$bkg), 1)
    expect_named(attributes(res)$bkg, c("A", "C", "G", "T"))
  }
  expect_error(run_mb(snps, motifs, threshold = 0.8, bkg = "bogus"), "must be")
})

test_that("parallel back-ends give the same results as serial", {
  snps <- example_variants()
  motifs <- hocomoco_core()
  serial <- example_run()
  key <- function(x) sort(paste(result_key(x), x$scoreRef, x$scoreAlt))

  skip_on_os("windows")
  ## regression: an empty result from one worker used to break the merge
  multicore <- run_mb(snps, motifs, filterp = TRUE, threshold = 1e-4,
                      BPPARAM = BiocParallel::MulticoreParam(2))
  expect_identical(key(multicore), key(serial))
})

test_that("SnowParam works and leaves a caller-started cluster running", {
  ## Snow workers load the installed package, so only test the installed copy
  skip_if(exists(".__DEVTOOLS__", envir = asNamespace("motifbreakR"), inherits = FALSE),
          "Snow workers would use the installed copy, not the one under test")
  snps <- example_variants()
  motifs <- hocomoco_core()
  serial <- example_run()
  key <- function(x) sort(paste(result_key(x), x$scoreRef, x$scoreAlt))

  p <- BiocParallel::SnowParam(2)
  BiocParallel::bpstart(p)
  on.exit(BiocParallel::bpstop(p), add = TRUE)
  snow <- run_mb(snps, motifs, filterp = TRUE, threshold = 1e-4, BPPARAM = p)
  expect_identical(key(snow), key(serial))
  expect_true(BiocParallel::bpisup(p))
})
