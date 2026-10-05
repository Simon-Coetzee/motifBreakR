## scoreMotif.R this module will score an arbitrary list of SNPs for motif
## disruption using a buffet of menu options for some sets of transcription
## factors

## optional: produce random SNP sets of equivalent size to the index SNP list to
## simulate background for enrichment calculations.

## define utility functions 'defaultOmega' (and related fxns) and maxPwm and
## minPwm

prepareVariants <- function(fsnplist, genome.bsgenome, max.pwm.width) {
  k <- as.integer(max.pwm.width); rm(max.pwm.width)
  ref_len <- nchar(fsnplist$REF)
  alt_len <- nchar(fsnplist$ALT)

  ## check that reference matches ref genome
  equals.ref <- getSeq(genome.bsgenome, fsnplist) == fsnplist$REF
  if (!all(equals.ref)) {
    stop(paste(names(fsnplist[!equals.ref]), "reference allele does not match value in reference genome. ",
               sep = " "))
  }

  ## centre every allele in a context of the same odd width; with an even
  ## maximum allele length, ref and alt contexts would otherwise differ in width
  max_len <- max(c(ref_len, alt_len))
  max_len <- max_len + (1L - max_len %% 2L)
  ref_pad <- as.integer((max_len - ref_len) %/% 2)
  alt_pad <- as.integer((max_len - alt_len) %/% 2)

  gr_ref <- fsnplist
  start(gr_ref) <- start(gr_ref) - k - ref_pad
  gr_ref <- resize(gr_ref, width = max_len + 2L * k, fix = "start")

  gr_alt <- fsnplist
  start(gr_alt) <- start(gr_alt) - k - alt_pad
  end(gr_alt) <- end(gr_alt) + k + alt_pad + (alt_len %% 2 == 0)

  snp.sequence.ref <- getSeq(genome.bsgenome, gr_ref)
  snp.sequence.alt <- getSeq(genome.bsgenome, gr_alt)

  at <- as(IRanges(start = k + alt_pad + 1L, width = ref_len), "IRangesList")
  snp.sequence.alt <- replaceAt(snp.sequence.alt, at, split(fsnplist$ALT, seq(fsnplist$ALT)))

  mcols(fsnplist)$contextRef <- snp.sequence.ref
  mcols(fsnplist)$contextCoordRef <- IRanges(start = k + ref_pad + 1L, end = k + ref_pad + ref_len)
  mcols(fsnplist)$contextAlt <- snp.sequence.alt
  mcols(fsnplist)$contextCoordAlt <- IRanges(start = k + alt_pad + 1L, end = k + alt_pad + alt_len)

  fsnplist$varType <- "Other"
  fsnplist$varType[width(fsnplist$REF) == 1L] <- "SNV"
  fsnplist$varType[width(fsnplist$REF) > width(fsnplist$ALT)] <- "Deletion"
  fsnplist$varType[width(fsnplist$REF) < width(fsnplist$ALT)] <- "Insertion"

  return(fsnplist)
}

varEff <- function(allelR, allelA) {
  score <- allelA - allelR
  effect <- cut(abs(score), breaks = c(-Inf, 0.4, 0.7, Inf), labels = c("neut", "weak", "strong"))
  names(effect) <- names(score)
  return(list(score = score, effect = effect))
}

reverseComplementMotif <- function(pwm) {
  rows <- rownames(pwm)
  cols <- colnames(pwm)
  Ns <- pwm["N", ,drop = FALSE]
  pwm <- pwm[4:1, length(cols):1, drop = FALSE]
  pwm <- rbind(pwm, Ns)
  rownames(pwm) <- rows
  colnames(pwm) <- cols
  return(pwm)
}


scoreSeqWindows <- function(ppm, seq, mask) {
  ppm.width <- ncol(ppm)
  seq.len <- nrow(seq)
  seq.w <- ncol(seq)
  recycle.w <-  1L + nrow(seq) - ppm.width
  ranges <- rep(seq.int(from = 1L, by = ppm.width,      to = recycle.w * ppm.width), each = ppm.width) +
                seq.int(from = 0L, by = ppm.width + 1L, length.out = ppm.width        )
  ranges <- rep.int(ranges, times = seq.w) + rep(seq.int(from = 0L, by = seq.len * ppm.width, length.out = seq.w), each = recycle.w * ppm.width)
  scores    <- colSums(matrix(t(                       ppm [seq, ])[ranges], nrow = ppm.width))
  scores_rc <- colSums(matrix(t(reverseComplementMotif(ppm)[seq, ])[ranges], nrow = ppm.width))
  res <- matrix(data = c(scores, scores_rc), ncol = recycle.w, byrow = TRUE)
  res_split <- rep(seq(seq.w), each = 2L)
  colnames(res) <- as.character(seq(ncol(res)))
  res <- res[res_split + c(0L,seq.w), ]
  rownames(res) <- rep.int(x = c(1L,2L), times = seq.w)
  res[seq(from = 1L, to = nrow(res), by = 2),][!mask] <- -Inf
  res[seq(from = 2L, to = nrow(res), by = 2),][!mask] <- -Inf
  res <- split.data.frame(res, res_split)
  names(res) <- colnames(seq)
  return(res)
}

maskWindows <- function(res, mw, nwindows, allele = "Ref") {
  rv <- length(res)

  col_idx <- matrix(seq(nwindows), nrow = rv, ncol = nwindows, byrow = TRUE)

  a_coord <- mcols(res)[[paste0("contextCoord", allele)]]
  a_start <- start(a_coord)
  a_end   <- end(a_coord)

  valid_mask <- (col_idx >= (a_start - mw + 1)) & (col_idx <= a_end)
  return(valid_mask)
}


maxThresholdWindows <- function(window.frame) {
  start.ind <- as.integer(colnames(window.frame)[1]) - 1L
  max.win <- arrayInd(which.max(window.frame), dim(window.frame))
  return(data.frame(window = as.integer(colnames(window.frame)[max.win[, 2] + start.ind]),
                    strand = c(1, 2)[max.win[, 1]]))
}

windowScore <- function(w.data) {
  return(w.data$score[w.data$strand, w.data$window])
}

passThresh <- function(ref.windows, alt.windows, thresh, filterp, pwmRanges) {
  if(filterp) {
    alt_pass <- vapply(alt.windows, function(x) any(x > thresh), logical(1))
    ref_pass <- vapply(ref.windows, function(x) any(x > thresh), logical(1))
  } else {
    alt_pass <- vapply(alt.windows, function(x) any((x - pwmRanges[1]) / (pwmRanges[2] - pwmRanges[1]) > thresh), logical(1))
    ref_pass <- vapply(ref.windows, function(x) any((x - pwmRanges[1]) / (pwmRanges[2] - pwmRanges[1]) > thresh), logical(1))
  }
  return(which(alt_pass | ref_pass))
}

#' Calculate Positions Relative to an Arbitrary Anchor (Vectorized)
#'
#' @noRd
#' @param result A GRanges object containing the results of a motif scan, with
#'   metadata columns for window index, context start and end, and motif strand
#'   for both alleles. The metadata columns should be named in the format
#'   "windowIdxRef", "contextStartRef", "contextEndRef", "motifStrandRef" for
#'   the reference allele and similarly for the alternate allele.
#' @param mw Integer (or vector). Motif width.
#' @param allele Character. "Ref" or "Alt" to specify which allele's context
#'   coordinates to use.
#' @param anchor_offset Integer. The offset from the context start to use as the
#'   anchor point (0-based).
#' @return A data.frame containing absolute and relative coordinates for all
#'   variants.
calculateAllPositions <- function(result, mw, allele, anchor_offset = 0L) {
  genomic_start        <- start(result)
  result               <- mcols(result)
  columns              <- c(toupper(allele), paste0(c("windowIdx", "motifStrand"), allele))
  context_start        <- 0L
  context_end          <- width(result[, columns[1]]) - 1L
  window_idx           <- result[, columns[2]]
  strand               <- result[, columns[3]]

  motif_context_start  <- window_idx
  motif_context_end    <- window_idx + mw - 1L

  anchor_idx           <- context_start + anchor_offset

  motif_rel_start      <- motif_context_start - anchor_idx
  motif_rel_end        <- motif_context_end - anchor_idx

  var_rel_start        <- context_start - anchor_idx
  var_rel_end          <- context_end - anchor_idx

  var_rel_motif_start  <- context_start - motif_context_start
  var_rel_motif_end    <- context_end - motif_context_start

  res <- DataFrame(
    motifCoordContext  = IRanges(start = motif_context_start, end = motif_context_end),
    motifCoordVariant  = IRanges(start = motif_rel_start, end = motif_rel_end),
    variantCoordOffset = IRanges(start = var_rel_start, end = var_rel_end),
    variantCoordMotif  = IRanges(start = var_rel_motif_start, end = var_rel_motif_end),
    motifStrand        = strand
  )

  colnames(res) <- paste0(colnames(res), allele)
  result[, colnames(res)] <- res

  return(result)
}

getMatchingSequence <- function(result, strongerIn) {
  res <- mcols(result)
  res$strongerIn <- strongerIn
  res$index <- seq(nrow(res))
  colquery <- c("motifCoordContext", "context")
  ref_res <- res[res$strongerIn == "Ref", c(paste0(colquery, "Ref"), "index")]
  colnames(ref_res) <- c(colquery, "index")
  alt_res <- res[res$strongerIn == "Alt", c(paste0(colquery, "Alt"), "index")]
  colnames(alt_res) <- c(colquery, "index")
  res <- rbind(ref_res, alt_res)
  res <- res[order(res$index), ]
  return(subseq(res$context, res$motifCoordContext))
}

#' @import methods
#' @import GenomicRanges
#' @import S4Vectors
#' @import BiocGenerics
#' @import IRanges
#' @importFrom Biostrings getSeq replaceLetterAt reverseComplement complement replaceAt matchPattern letterFrequency subseq
#' @importFrom pwalign pairwiseAlignment insertion deletion
#' @importFrom matrixStats colAlls
#' @importFrom TFMPvalue TFMpv2sc
#' @importFrom stringr str_locate_all str_sub
scoreSnpList <- function(pwmList, fsnplist, method = "default", bkg = NULL,
                         show.neutral = FALSE, verbose = FALSE,
                         genome.bsgenome = NULL, filterp = TRUE) {

  threshold <- pwmList$pwmThreshold
  pwmList.pc <- pwmList$pwmListPseudoCount
  pwmRanges <- pwmList$pwmRange
  pwmList <- pwmList$pwmList
  pwmConsensus <- DataFrame(
    consensus = DNAStringSet(lapply(pwmList.pc, function(pwm) {
      DNAString(paste0(rownames(pwm)[apply(pwm , 2, function(x) {which.max(x)})], collapse = "")) })),
    revComConsensus = DNAStringSet(lapply(pwmList.pc, function(pwm) {
      DNAString(paste0(rownames(pwm)[apply(reverseComplementMotif(pwm), 2, function(x) {which.max(x)})], collapse = "")) })))

  k <- max(vapply(pwmList, ncol, integer(1)))


  snp.sequence.alt <- fsnplist$contextAlt
  snp.sequence.ref <- fsnplist$contextRef
  fsnplist$index <- seq(fsnplist)
  results <- rep(fsnplist, times = length(pwmList))

  mcols(results) <- c(mcols(results), rep(DataFrame(windowIdxRef = integer(length(pwmList)),
                                                    motifStrandRef = "*",
                                                    windowIdxAlt = integer(length(pwmList)),
                                                    motifStrandAlt = "*",
                                                    motifID = mcols(pwmList)$providerId,
                                                    geneSymbol = mcols(pwmList)$geneSymbol,
                                                    dataSource = mcols(pwmList)$dataSource,
                                                    providerName = mcols(pwmList)$providerName,
                                                    providerId = mcols(pwmList)$providerId,
                                                    seqMatch = character(length(pwmList)),
                                                    pwmConsensus = character(length(pwmList)),
                                                    pctRef = numeric(length(pwmList)),
                                                    pctAlt = numeric(length(pwmList)),
                                                    scoreRef = numeric(length(pwmList)),
                                                    scoreAlt = numeric(length(pwmList)),
                                                    pValueRef = NA_real_,
                                                    pValueAlt = NA_real_,
                                                    strongerIn = character(length(pwmList)),
                                                    alleleDiff = numeric(length(pwmList)),
                                                    alleleEffectSize = numeric(length(pwmList)),
                                                    effect = numeric(length(pwmList))),
                                          each = length(fsnplist)))
  if (!filterp) {
    mcols(results)[, c("pValueRef", "pValueAlt")] <- NULL
  }

  ref.len <- width(fsnplist$REF)
  alt.len <- width(fsnplist$ALT)

  snp.sequence.ref.pwm <- t(as.matrix(snp.sequence.ref))
  snp.sequence.alt.pwm <- t(as.matrix(snp.sequence.alt))
  resultSet <- GRanges()
  for (pwm.i in seq_along(pwmList)) {
    pwm.basic <- pwmList[[pwm.i]]
    pwm <- pwmList.pc[[pwm.i]]
    thresh <- threshold[[pwm.i]]
    pwm.len <- ncol(pwm)
    n.window <- width(snp.sequence.ref[1]) - pwm.len + 1L

    ref_mask <- maskWindows(fsnplist, mw = pwm.len, nwindows = n.window, allele = "Ref")
    alt_mask <- maskWindows(fsnplist, mw = pwm.len, nwindows = n.window, allele = "Alt")

    ref.windows <- scoreSeqWindows(ppm = pwm, seq = snp.sequence.ref.pwm, mask = ref_mask)
    alt.windows <- scoreSeqWindows(ppm = pwm, seq = snp.sequence.alt.pwm, mask = alt_mask)

    pass_effect <- passThresh(ref.windows, alt.windows, thresh, filterp, pwmRanges[[pwm.i]])
    if (any(!is.na(pass_effect))) {
      ref.windows <- ref.windows[pass_effect]
      alt.windows <- alt.windows[pass_effect]
      hit.alt <- do.call(rbind, lapply(alt.windows, maxThresholdWindows))
      hit.ref <- do.call(rbind, lapply(ref.windows, maxThresholdWindows))
      ref.len <- ref.len[pass_effect]
      alt.len <- alt.len[pass_effect]

      allelR <- mapply(function(windows, strand, window) {
        windows[strand, window]
      }, windows = ref.windows, strand = hit.ref$strand, window = hit.ref$window)
      allelA <- mapply(function(windows, strand, window) {
        windows[strand, window]
      }, windows = alt.windows, strand = hit.alt$strand, window = hit.alt$window)

      scorediff <- varEff(allelR, allelA)
      effect <- scorediff$effect
      score <- scorediff$score

      if (show.neutral) {
        keep_index <- seq_along(pass_effect)
      } else {
        keep_index <- which(effect != "neut")
      }
      motif_keep <- pass_effect[keep_index]

      if(all(is.na(motif_keep)) | length(keep_index) < 1) next()

      ## results holds one block of length(fsnplist) rows per PWM; select by
      ## position since providerId is not unique across motif collections
      motif_result <- results[(pwm.i - 1L) * length(fsnplist) + motif_keep]

      motif_result$windowIdxRef <- hit.ref[keep_index,]$window
      motif_result$motifStrandRef <- factor(hit.ref[keep_index,]$strand, levels = c(1,2), labels = c("+", "-"))
      motif_result$windowIdxAlt <- hit.alt[keep_index,]$window
      motif_result$motifStrandAlt <- factor(hit.alt[keep_index,]$strand, levels = c(1,2), labels = c("+", "-"))

      mcols(motif_result) <- calculateAllPositions(result = motif_result,
                                                   mw = pwm.len,
                                                   allele = "Ref",
                                                   anchor_offset = 0L)
      mcols(motif_result) <- calculateAllPositions(result = motif_result,
                                                   mw = pwm.len,
                                                   allele = "Alt",
                                                   anchor_offset = 0L)

      rev_com_motif <- ifelse(score[keep_index] > 0,
                              motif_result$motifStrandAlt,
                              motif_result$motifStrandRef)
      mDF <- DataFrame(pctRef = (allelR[keep_index] - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1]),
                       pctAlt = (allelA[keep_index] - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1]),
                       scoreRef = allelR[keep_index],
                       scoreAlt = allelA[keep_index],
                       windowIdxRef = motif_result$windowIdxRef - start(motif_result$contextCoordRef),
                       windowIdxAlt = motif_result$windowIdxAlt - start(motif_result$contextCoordAlt),
                       alleleDiff = score[keep_index],
                       strongerIn = ifelse(score[keep_index] > 0, "Alt", "Ref"),
                       pwmConsensus = c(pwmConsensus$consensus[pwm.i],
                                        pwmConsensus$revComConsensus[pwm.i])[as.integer(rev_com_motif)],
                       seqMatch = getMatchingSequence(motif_result, ifelse(score[keep_index] > 0, "Alt", "Ref")),
                       alleleEffectSize = score[keep_index]/pwmRanges[[pwm.i]][[2]],
                       effect = effect[keep_index])
      mcols(motif_result)[,colnames(mDF)] <- mDF
      mcols(motif_result)[, c("motifCoordContextRef", "motifCoordVariantRef",
                              "variantCoordOffsetRef", "variantCoordMotifRef",
                              "contextRef", "contextCoordRef",
                              "motifCoordContextAlt", "motifCoordVariantAlt",
                              "variantCoordOffsetAlt", "variantCoordMotifAlt",
                              "contextAlt", "contextCoordAlt",
                              "index")] <- NULL
      resultSet <- c(resultSet, motif_result)
    }
  }
  if (length(resultSet) < 1) {
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist),
                    "with 0 potentially disruptive matches to", length(unique(resultSet$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
    return(NULL)
  } else {
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist),
                    "with", length(resultSet), "potentially disruptive matches to", length(unique(resultSet$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
    return(resultSet)
  }
}

#' @importFrom matrixStats colRanges
#' @importFrom stringr str_pad
#' @importFrom TFMPvalue TFMsc2pv
#' @importFrom matrixStats colMaxs colMins
preparePWM <- function(pwmList,
                       filterp,
                       bkg,
                       scoreThresh,
                       method = "default") {

  bkg <- bkg[c('A', 'C', 'G', 'T')]
  scounts <- as.integer(mcols(pwmList)$sequenceCount)
  mCGMAT <- pwmList@manuallyCuratedGeneMotifAssociationTable
  scounts[is.na(scounts)] <- 20L
  pwmList.pc <- Map(function(pwm, scount) {
    pwm <- (pwm * scount + bkg)/(scount + 1)
  }, pwmList, scounts)
  if (method == "ic") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm, b=bkg) {
      omegaic <- colSums(pwm * log2(pwm/b))
    })
  }
  if (method == "default") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm) {
      omegadefault <- colMaxs(pwm) - colMins(pwm)
    })
  }
  if (method == "log") {
    pwmList.pc <- lapply(pwmList.pc, function(pwm, b) {
      pwm <- log2(pwm) - log2(b)
    }, b = bkg)
    pwmOmegas <- 1
  }
  if (method == "notrans") {
    pwmOmegas <- 1
  }
  pwmList.pc <- Map(function(pwm, omega) {
    if (length(omega) == 1 && omega == 1) {
      return(pwm)
    } else {
      omegamatrix <- matrix(rep(omega, 4), nrow = 4, byrow = TRUE)
      pwm <- pwm * omegamatrix
    }
  }, pwmList.pc, pwmOmegas)
  pwmRanges <- Map(function(pwm, omega) {
    x <- colSums(colRanges(pwm))
    return(x)
  }, pwmList.pc, pwmOmegas)
  if (filterp) {
    pwmList.pc2 <- lapply(pwmList.pc, round, digits = 2)
    pwmThresh <- lapply(pwmList.pc2, TFMpv2sc, pvalue = scoreThresh, bg = bkg, type = "PWM")
    pwmThresh <- Map("+", pwmThresh, -0.02)
  } else {
    pwmThresh <- rep.int(scoreThresh, times = length(pwmRanges))
  }
  pwmList@listData <- lapply(pwmList, function(pwm) {
    pwm <- pwm[c("A", "C", "G", "T"), ]
    pwm <- rbind(pwm, N = 0)
    colnames(pwm) <- as.character(1:ncol(pwm))
    return(pwm) })
  pwmList.pc <- lapply(pwmList.pc, function(pwm) {
    pwm <- pwm[c("A", "C", "G", "T"), ]
    pwm <- rbind(pwm, N = 0)
    colnames(pwm) <- as.character(1:ncol(pwm))
    return(pwm) })
  return(list(pwmList = pwmList,
              pwmListPseudoCount = pwmList.pc,
              pwmRange = pwmRanges,
              pwmThreshold = pwmThresh))
}

get_background <- function(bkg, snpList, genome.bsgenome, pwmList) {
  if(is.character(bkg) && bkg == "genome") {
    bg <- colSums(letterFrequency(getSeq(genome.bsgenome), letters = c("A", "C", "G", "T")))
  } else if (is.character(bkg) && bkg == "aggregate") {
    bg <- colSums(letterFrequency(snpList$contextRef, letters = c("A", "C", "G", "T")))
  } else if (is.numeric(bkg) && length(bkg) == 4) {
    bg <- bkg
  } else if (is.character(bkg) && bkg == "pwm") {
    bg <- rowSums(sapply(pwmList, rowSums))[c("A", "C", "G", "T")]
  } else {
    stop("bg must be 'genome', 'aggregate', 'pwm', or a numeric vector of length 4")
  }
  bg <- bg / sum(bg)
  return(bg)
}

#' Predict The Disruptiveness Of Genetic Variants On Transcription Factor
#' Binding Sites
#'
#' @param snpList The output of \code{\link{snps.from.rsid}} or
#'   \code{\link{snps.from.file}}; may contain SNVs and indels
#' @param pwmList An object of class \code{MotifList} containing the motifs that
#'   you wish to interrogate
#' @param threshold Numeric; with \code{filterp = TRUE}, the maximum p-value for
#'   a match to be reported; otherwise the minimum score, as a fraction (0-1) of
#'   the motif's scoring range (see \code{pctRef} and \code{pctAlt}).
#' @param method Character; one of \code{default}, \code{log}, \code{ic}, or
#'   \code{notrans}; see Details.
#' @param bkg the background nucleotide frequencies; either a numeric vector of
#'   length 4 named \code{A}, \code{C}, \code{G} and \code{T}, or one of:
#'   \describe{
#'     \item{\code{"genome"}}{frequencies across the whole reference genome. This
#'       reads the entire genome sequence into memory.}
#'     \item{\code{"aggregate"}}{frequencies across the sequence windows
#'       surrounding the input variants.}
#'     \item{\code{"pwm"}}{the average nucleotide composition of the motifs in
#'       \code{pwmList}.}
#'   }
#'   The background is used for the motif pseudocounts, by \code{method = "log"}
#'   and \code{method = "ic"}, and for p-values. It is stored in
#'   \code{attributes(results)$bkg} and reused by \code{\link{calculatePvalue}}.
#' @param filterp Logical; filter by p-value instead of by pct score.
#' @param show.neutral Logical; include neutral changes in the output
#' @param verbose Logical; show progress messages (from the workers only when
#'   running serially)
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam-class}} object
#'   controlling parallel evaluation; work is split across workers by motif.
#'   Try \code{BiocParallel::registered()} to see what is available; for example
#'   \code{BiocParallel::SerialParam()} gives serial evaluation, and
#'   \code{BiocParallel::SnowParam()} parallel evaluation on Windows.
#' @seealso See \code{\link{snps.from.rsid}} and \code{\link{snps.from.file}} for
#'   information about how to generate the input to this function and
#'   \code{\link{plotMB}} for information on how to visualize its output
#' @details \pkg{motifbreakR} works with position probability matrices (PPM). PPM
#' are derived as the fractional occurrence of nucleotides A,C,G, and T at
#' each position of a position frequency matrix (PFM). PFM are simply the
#' tally of each nucleotide at each position across a set of aligned
#' sequences. With a PPM, one can generate probabilities based on the
#' genome, or more practically, create any number of position specific
#' scoring matrices (PSSM) based on the principle that the PPM contains
#' information about the likelihood of observing a particular nucleotide at
#' a particular position of a true transcription factor binding site. What
#' follows is a discussion of the different algorithms that may be
#' employed in calls to the \pkg{motifbreakR} function via the \code{method}
#' argument.
#'
#' Before scoring, a pseudocount is added to each PPM so that logarithms are
#' defined: \code{ppm <- (ppm * sequenceCount + bkg) / (sequenceCount + 1)},
#' where \code{sequenceCount} is taken from the motif metadata, or 20 when it is
#' \code{NA}.
#'
#' Suppose we have a frequency matrix \eqn{M} of width \eqn{n} (\emph{i.e.} a
#' PPM as described above). Furthermore, we have a sequence \eqn{s} also of
#' length \eqn{n}, such that
#' \eqn{s_{i} \in \{ A,T,C,G \}, i = 1,\ldots n}{s_i in {A,T,G,C}, i = 1 \ldots n}.
#' Each column of
#' \eqn{M} contains the frequencies of each letter in each position.
#'
#' Commonly in the literature sequences are scored as the sum of log
#' probabilities:
#'
#' \strong{Equation 1}
#'
#' \deqn{F( s,M ) = \sum_{i = 1}^{n}{\log_2( \frac{M_{s_{i},i}}{b_{s_{i}}} )}}{
#' F( s,M ) = \sum_(i = 1)^n log2 ((M_s_i,_i)/b_s_i)}
#'
#' where \eqn{b_{s_{i}}}{b_s_i} is the background frequency of letter \eqn{s_{i}}{s_i} in
#' the genome of interest. This method can be specified by the user as
#' \code{method='log'}.
#'
#' As an alternative to this method, we introduced a scoring method to
#' directly weight the score by the importance of the position within the
#' match sequence. This method of weighting is accessed by specifying
#' \code{method='default'} or \code{method='ic'}. A general representation
#' of this scoring method is given by:
#'
#' \strong{Equation 2}
#'
#' \deqn{F( s,M ) = p_{s} \cdot \omega_{M}}{F( s,M ) = p_s . \omega_M}
#'
#' where \eqn{p_{s}}{p_s} is the scoring vector derived from sequence \eqn{s} and matrix
#' \eqn{M}, and \eqn{w_{M}}{w_M} is a weight vector derived from \eqn{M}. First, we
#' compute the scoring vector of position scores \eqn{p}
#'
#' \strong{Equation 3}
#'
#' \deqn{p_{s} = ( M_{s_{i},i} ) \textrm{\ \ \ for\ \ \ } i = 1,\ldots n}{
#' p_s = ( M_s_i,_i ) for i = 1 \ldots n}
#'
#' and second, for each \eqn{M} a constant vector of weights
#' \eqn{\omega_{M} = ( \omega_{1},\omega_{2},\ldots,\omega_{n} )}{\omega_M = ( \omega_1, \omega_2, \ldots, \omega_n)}.
#'
#' There are two methods for producing \eqn{\omega_{M}}{\omega_M}. The first, which we
#' call weighted sum (\code{method='default'}), is the difference of the maximum
#' and minimum values for each column of \eqn{M}:
#'
#' \strong{Equation 4.1}
#'
#' \deqn{\omega_{i} = \max \{ M_{i} \} - \min \{ M_{i} \}\textrm{\ \ \ \ where\ \ \ \ \ \ }i = 1,\ldots n}{
#' \omega_i = max{M_i} - min{M_i} where i = 1 \ldots n}
#'
#' The second variation of this theme is to weight by relative entropy.
#' Thus the relative entropy weight for each column \eqn{i} of the matrix is
#' given by:
#'
#' \strong{Equation 4.2}
#'
#' \deqn{\omega_{i} = \sum_{j \in \{ A,C,G,T \}}^{}{M_{j,i}\log_2( \frac{M_{j,i}}{b_{j}} )}\textrm{\ \ \ \ \ where\ \ \ \ \ }i = 1,\ldots n}{
#' \omega_i = \sum_{j in {A,C,G,T}} {M_(j,i)} log2(M_(j,i)/b_j) where i = 1 \ldots n}
#'
#' where \eqn{b_{j}}{b_j} is again the background frequency of the letter \eqn{j}.
#'
#' Thus, there are 3 scoring algorithms to apply via the \code{method}
#' argument. The first is the standard summation of log probabilities
#' (\code{method='log'}). The second and third are the weighted sum and
#' information content methods (\code{method='default'} and \code{method='ic'}) specified by
#' equations 4.1 and 4.2, respectively. Additionally \code{method='notrans'}
#' scores the sum of the untransformed probabilities, \eqn{\sum p_{s}}{sum(p_s)}. \pkg{motifbreakR} assumes a
#' uniform background nucleotide distribution (\eqn{b}) in equations 1 and
#' 4.2 unless otherwise specified by the user. Since we are primarily
#' interested in the difference between alleles, background frequency is
#' not a major factor, although it can change the results. Additionally,
#' inclusion of background frequency introduces potential bias when
#' collections of motifs are employed, since motifs are themselves
#' unbalanced with respect to nucleotide composition. With these cautions
#' in mind, users may override the uniform distribution if so desired. For
#' all three methods, \pkg{motifbreakR} scores and reports the reference
#' and alternate alleles of the sequence
#' (\eqn{F( s_{\textsc{ref}},M )}{F( s_ref,M )} and
#' \eqn{F( s_{\textsc{alt}},M )}{F( s_alt,M )}), and provides the matrix scores
#' \eqn{p_{s_{\textsc{ref}}}}{p_s_ref} and \eqn{p_{s_{\textsc{alt}}}}{p_s_alt} of the SNP (or
#' variant). The scores are scaled as a fraction of scoring range 0-1 of
#' the motif matrix, \eqn{M}. If either of
#' \eqn{F( s_{\textsc{ref}},M )}{F( s_ref,M )} and
#' \eqn{F( s_{\textsc{alt}},M )}{F( s_alt,M )} is greater than a user-specified
#' threshold (default value of 0.85) the SNP is reported. The effect of a
#' variant is classified from the absolute difference in score between the
#' alleles (\code{alleleDiff}): \code{"strong"} above 0.7, \code{"weak"} from
#' 0.4 to 0.7 and \code{"neut"} below 0.4. By default \pkg{motifbreakR} does not
#' report neutral effects; set \code{show.neutral = TRUE} to include them.
#'
#' Additionally, now, with the use of \code{\link[TFMPvalue]{TFMPvalue-package}}, we may filter by p-value of the match.
#' This is unfortunately a two step process. First, by invoking \code{filterp=TRUE} and setting a threshold at
#' a desired p-value e.g 1e-4, we perform a rough filter on the results by rounding all values in the PWM to two
#' decimal places, and calculating a scoring threshold based upon that. The second step is to use the function \code{\link{calculatePvalue}()}
#' on a selection of results which will change the \code{pValueRef} and \code{pValueAlt} columns in the output from \code{NA} to the p-value
#' calculated by \code{\link[TFMPvalue]{TFMsc2pv}}.  This can be (although not always) a very memory and time intensive process if the algorithm doesn't converge rapidly.
#'
#' @return a GRanges object containing:
#'  \item{SNP_id}{the identifier of the variant}
#'  \item{REF}{the reference allele for the variant}
#'  \item{ALT}{the alternate allele for the variant}
#'  \item{varType}{one of \code{SNV}, \code{Insertion}, \code{Deletion}, or \code{Other}}
#'  \item{windowIdxRef, windowIdxAlt}{the start of the best motif match on the
#'  reference and alternate allele, relative to the first base of the variant}
#'  \item{motifStrandRef, motifStrandAlt}{the strand (\code{+} or \code{-}) of the
#'  best motif match on the reference and alternate allele}
#'  \item{geneSymbol}{the geneSymbol corresponding to the TF of the TF binding motif}
#'  \item{dataSource}{the source of the TF binding motif}
#'  \item{motifID, providerName, providerId}{the name and id provided by the source}
#'  \item{seqMatch}{the sequence matched by the TF binding motif on the allele
#'  with the stronger match (see \code{strongerIn})}
#'  \item{pwmConsensus}{the consensus sequence of the TF binding motif, on the
#'  strand of the stronger match}
#'  \item{pctRef}{The score as determined by the scoring method, when the sequence contains the reference variant allele, normalized to a scale from 0 - 1. If \code{filterp = FALSE},
#'  this is the value that is thresholded.}
#'  \item{pctAlt}{The score as determined by the scoring method, when the sequence contains the alternate variant allele, normalized to a scale from 0 - 1. If \code{filterp = FALSE},
#'  this is the value that is thresholded.}
#'  \item{scoreRef}{The score as determined by the scoring method, when the sequence contains the reference variant allele}
#'  \item{scoreAlt}{The score as determined by the scoring method, when the sequence contains the alternate variant allele}
#'  \item{pValueRef}{p-value for the match for the pctRef score, initially set to \code{NA}; only present when \code{filterp = TRUE}. see \code{\link{calculatePvalue}} for more information}
#'  \item{pValueAlt}{p-value for the match for the pctAlt score, initially set to \code{NA}; only present when \code{filterp = TRUE}. see \code{\link{calculatePvalue}} for more information}
#'  \item{strongerIn}{\code{Ref} or \code{Alt}, the allele with the stronger motif match}
#'  \item{alleleDiff}{The difference between the score on the reference allele and the score on the alternate allele}
#'  \item{alleleEffectSize}{The ratio of the \code{alleleDiff} and the maximal score of a sequence under the PWM}
#'  \item{effect}{one of \code{"strong"}, \code{"weak"}, or \code{"neut"}
#'  (neutral, only with \code{show.neutral = TRUE}) indicating the strength of the
#'  effect; see Details.}
#'  each SNP in this object may be plotted with \code{\link{plotMB}}
#' @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  # prepare variants
#'  load(system.file("extdata",
#'                   "pca.enhancer.snps.rda",
#'                   package = "motifbreakR")) # loads snps.mb
#'  pca.enhancer.snps <- sample(snps.mb, 20)
#'  # Get motifs to interrogate
#'  data(hocomoco)
#'  motifs <- sample(hocomoco, 50)
#'  # run motifbreakR
#'  results <- motifbreakR(pca.enhancer.snps,
#'                         motifs, threshold = 0.85,
#'                         method = "ic",
#'                         BPPARAM=BiocParallel::SerialParam())
#' @import BiocParallel
#' @importClassesFrom MotifDb MotifList
#' @importFrom BiocParallel bplapply bpnworkers
#' @importFrom stringr str_length str_trim
#' @export
motifbreakR <- function(snpList, pwmList, threshold = 0.85, filterp = FALSE,
                        method = "default", show.neutral = FALSE, verbose = FALSE,
                        bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                        BPPARAM = bpparam()) {

  if (.Platform$OS.type == "windows" && inherits(BPPARAM, "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n",
            "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpnworkers(BPPARAM)
  num.snps <- length(snpList)
  if (num.snps < cores) {
    cores <- num.snps
  }

  genome.package <- attributes(snpList)$genome.package
  if (requireNamespace(genome.package, quietly = TRUE)) {
    genome.bsgenome <- getExportedValue(genome.package, genome.package)
  } else {
    stop(paste0(genome.package, " is the genome selected for this snp list and \n",
                "  is not present on your environment. Please load it and try again."))
  }

  k <- max(vapply(pwmList, ncol, integer(1)))

  snpList <- prepareVariants(fsnplist = snpList,
                             genome.bsgenome = genome.bsgenome,
                             max.pwm.width = k)

  bg <- get_background(bkg, snpList, genome.bsgenome, pwmList)

  split_cores <- split(names(pwmList), rep_len(sequence(cores), length(pwmList)))
  pwms <- lapply(split_cores, function(pwm_i) {
    preparePWM(pwmList = pwmList[pwm_i], filterp = filterp,
               scoreThresh = threshold, bkg = bg,
               method = method)
  })

  x <- bplapply(pwms, scoreSnpList,
                fsnplist = snpList,
                method = method,
                bkg = bg,
                show.neutral = show.neutral,
                verbose = ifelse(cores == 1, verbose, FALSE),
                genome.bsgenome = genome.bsgenome,
                filterp = filterp,
                BPPARAM = BPPARAM)

  pwms <- preparePWM(pwmList = pwmList, filterp = filterp,
                     scoreThresh = threshold, bkg = bg,
                     method = method)

  pwmList <- pwms$pwmList
  pwmList@listData <- lapply(pwms$pwmList, function(pwm) { pwm <- pwm[c("A", "C", "G", "T"), ]; return(pwm) })
  pwmList.pc <- lapply(pwms$pwmListPseudoCount, function(pwm) { pwm <- pwm[c("A", "C", "G", "T"), ]; return(pwm) })

  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) > 0) {
    x <- unlist(GRangesList(unname(x)))
    x <- x[order(match(x$SNP_id, names(snpList)), x$geneSymbol), ]
    attributes(x)$genome.package <- genome.package
    attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% unique(x$providerId) &
                                      mcols(pwmList)$providerName %in% unique(x$providerName), ]
    attributes(x)$scoremotifs <- pwmList.pc[names(attributes(x)$motifs)]
    attributes(x)$bkg <- bg
  } else {
    warning("No SNP/Motif Interactions reached threshold")
    x <- NULL
  }
  if (verbose) {
    if (is.null(x)) {
      message(paste("reached end of SNPs list length =", num.snps, "with 0 potentially disruptive matches to",
                    length(unique(x$geneSymbol)), "of", length(pwmList), "motifs."))

    } else {
      message(paste("reached end of SNPs list length =", num.snps, "with",
                    length(x), "potentially disruptive matches to", length(unique(x$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
  }
  return(x)
}


#' Calculate the significance of the matches for the reference and alternate alleles for their PWM
#'
#' @param results The output of \code{motifbreakR} that was run with \code{filterp=TRUE}
#' @param granularity Numeric; the granularity to which to round the PWM,
#'  larger values compromise full accuracy for speed of calculation. A value of
#'  \code{NULL} does no rounding.
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam-class}} object
#'   controlling parallel evaluation. Try \code{BiocParallel::registered()} to
#'   see what is available; the default \code{BiocParallel::SerialParam()} gives
#'   serial evaluation.
#' @return a GRanges object. The same GRanges object that was input as \code{results}, but with
#'  \code{pValueRef} and \code{pValueAlt} columns in the output modified from \code{NA} to the p-value
#'  calculated by \code{\link[TFMPvalue]{TFMsc2pv}}. Additionally a \code{pValueEffect} column that indicates "strong"
#'  when the two p-values differ by more than an order of magnitude, otherwise "weak".
#' @seealso See \code{\link[TFMPvalue]{TFMsc2pv}} from the \pkg{TFMPvalue} package for
#'   information about how the p-values are calculated.
#' @details This function is intended to be used on a selection of results produced by \code{\link{motifbreakR}}, and
#' this can be (although not always) a very memory and time intensive process if the algorithm doesn't converge rapidly.
#' The nucleotide background used is the one stored by \code{motifbreakR} in
#' \code{attributes(results)$bkg}.
#' @source Hélène Touzet and Jean-Stéphane Varré (2007) Efficient and accurate P-value computation for Position Weight Matrices.
#'  Algorithms for Molecular Biology, \bold{2: 15}.
#' @examples
#' data(example.results)
#' rs1006140 <- example.results[example.results$SNP_id %in% "rs1006140"]
#' # low granularity for speed; 1e-6 or 1e-7 recommended for accuracy
#' rs1006140 <- calculatePvalue(rs1006140, BPPARAM=BiocParallel::SerialParam(), granularity = 1e-4)
#'
# #' @importFrom qvalue qvalue
#'
#' @export
calculatePvalue <- function(results,
                            granularity = NULL,
                            BPPARAM = BiocParallel::SerialParam()) {

  if (.Platform$OS.type == "windows" && inherits(BPPARAM, "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n",
                   "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpnworkers(BPPARAM)
  num.res <- length(results)
  if (num.res < cores) {
    cores <- num.res
  }
  if(!("pValueRef" %in% names(mcols(results)))) {
    stop('incorrect results format; please rerun analysis with filterp=TRUE')
  } else {
    pwmListmeta <- mcols(attributes(results)$motifs, use.names=TRUE)
    pwmList <- attributes(results)$scoremotifs
    background <- attributes(results)$bkg
    if(!is.null(granularity)) {
      pwmList <- lapply(pwmList, function(x, g) {x <- floor(x/g)*g; return(x)}, g = granularity)
    }
    results_sp <- split(results, seq_along(results))
    pvalues <- bplapply(results_sp, function(i, pwmList, pwmListmeta, bkg) {
      result <- i
      pwm.id <- result$providerId
      pwm.name.f <- result$providerName
      pwmmeta <- pwmListmeta[pwmListmeta$providerId == pwm.id & pwmListmeta$providerName == pwm.name.f, ]
      pwm <- pwmList[[rownames(pwmmeta)[1]]]
      ref <- TFMsc2pv(pwm, mcols(result)[["scoreRef"]], bg = bkg, type="PWM")
      alt <- TFMsc2pv(pwm, mcols(result)[["scoreAlt"]], bg = bkg, type="PWM")
      gc()
      return(data.frame(ref=ref, alt=alt))
    }, pwmList=pwmList, pwmListmeta=pwmListmeta, bkg = background, BPPARAM = BPPARAM)
    pvalues.df <- base::do.call("rbind", c(pvalues, make.row.names = FALSE))
    results$pValueRef <- pvalues.df[, "ref"]
    results$pValueAlt <- pvalues.df[, "alt"]
    pscore <- with(results, ifelse(pValueRef < pValueAlt, pValueAlt/pValueRef, pValueRef/pValueAlt))
    results$pValueEffect <- ifelse(pscore > 10, "strong", "weak")

    return(results)
  }
}

addPWM.stack <- function(identifier, index, GdObject, pwm_stack, ...) {
  plotMotifLogoStack.3(pwm_stack)
}

selcor <- function(identifier, index, GdObject, ... ) {
  if (identical(index, 1L)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

selall <- function(identifier, GdObject, ... ) {
    return(TRUE)
}

#' @importFrom grid grid.newpage pushViewport viewport popViewport
plotMotifLogoStack.3 <- function(pfms, ...) {
  n <- length(pfms)
  lapply(pfms, function(.ele) {
    if (!is(.ele, 'pfm'))
      stop("pfms must be a list of class pfm")
  })
  assign("tmp_motifStack_symbolsCache", list(), pos = ".GlobalEnv")
  ht <- 1/n
  y0 <- 0.5 * ht
  for (i in rev(seq.int(n))) {
    pushViewport(viewport(y = y0, height = ht))
    suppressWarnings(plotMotifLogo(pfms[[i]], motifName = pfms[[i]]@name, ncex = 1,
                                   p = pfms[[i]]@background, colset = pfms[[i]]@color,
                                   xlab = NA, newpage = FALSE, margins = c(1.5, 4.1,
                                                                           1.1, 0.1), ...))
    popViewport()
    y0 <- y0 + ht
  }
  rm(list = "tmp_motifStack_symbolsCache", pos = ".GlobalEnv")
  return()
}

getAlleleSpecific <- function(result, col, coord = "start", offset = 0L) {
  if (!(coord %in% c("start", "end", "seq"))) {
    stop("coord must be one of 'start', 'end', or 'seq'")
  }

  condition <- result$alleleDiff > 0

  mcols(result) <- calculateAllPositions(result = result, mw = width(result$pwmConsensus), allele = "Ref", anchor_offset = offset)
  mcols(result) <- calculateAllPositions(result = result, mw = width(result$pwmConsensus), allele = "Alt", anchor_offset = offset)

  alt_col <- paste0(col, "Alt")
  ref_col <- paste0(col, "Ref")

  if (coord == "seq" ) {
    fun <- c
    out <- rep.int(IRanges(start = 0, width = 1), length(result))
  } else {
    fun <- getExportedValue("BiocGenerics", coord)
    out <- numeric(length(result))
  }

  out[condition] <- fun(mcols(result)[condition, alt_col])
  out[!condition] <- fun(mcols(result)[!condition, ref_col])

  return(out)
}

getGaps <- function(result) {
  length_diffs <- lengths(result$ALT) - lengths(result$REF)
  add_gaps <- (-length_diffs * (result$alleleDiff > 0 & result$varType == "Deletion")) +
    (length_diffs * (result$alleleDiff < 0 & result$varType == "Insertion"))
  return(add_gaps)
}

#' @importFrom stringr str_remove
#' @importFrom motifStack addBlank
DNAmotifAlignment.2snp <- function(pwms, result) {
  length_diffs <- lengths(result$ALT) - lengths(result$REF)
  add_gaps <- (-length_diffs * (result$alleleDiff > 0 & result$varType == "Deletion")) +
    (length_diffs * (result$alleleDiff < 0 & result$varType == "Insertion"))

  froms <- getAlleleSpecific(result, "motifCoordVariant", coord = "start")
  tos <- getAlleleSpecific(result, "motifCoordVariant", coord = "end")
  from = min(froms)
  to = max(tos + add_gaps)

  for (pwm.i in seq_along(pwms)) {
    pwm.name <- pwms[[pwm.i]]@name
    pwm.name <- str_remove(pwm.name, pattern = "-:rc$|-:r$")

    pwm.info <- attributes(result)$motifs
    pwm.id <- mcols(pwm.info[pwm.name, ])$providerId
    pwm.name <- mcols(pwm.info[pwm.name, ])$providerName
    pwm.w <- ncol(pwms[[pwm.i]]@mat)
    mresult <- result[result$providerId == pwm.id & result$providerName == pwm.name, ]

    length_diff <- lengths(mresult$ALT) - lengths(mresult$REF)
    add_gap <- (-length_diff * (mresult$alleleDiff > 0 & mresult$varType == "Deletion")) +
      (length_diff * (mresult$alleleDiff < 0 & mresult$varType == "Insertion"))

    mstart <- start(intersect(IRanges(0, pwm.w), getAlleleSpecific(mresult, "variantCoordMotif", coord = "seq"))) + 1L
    new.mat <- cbind(pwms[[pwm.i]]@mat[, seq(mstart), drop = FALSE],
                     matrix(0.25, nrow = 4, ncol = add_gap),
                     pwms[[pwm.i]]@mat[, seq(from = mstart + 1L, to = pwm.w), drop = FALSE])
    pwms[[pwm.i]]@mat <- new.mat

    start.offset <- getAlleleSpecific(mresult, "motifCoordVariant", coord = "start") - from
    end.offset <- to - (getAlleleSpecific(mresult, "motifCoordVariant", coord = "end") + add_gap)
    pwms[[pwm.i]] <- addBlank(x = pwms[[pwm.i]], n = start.offset, b = FALSE)
    pwms[[pwm.i]] <- addBlank(x = pwms[[pwm.i]], n = end.offset, b = TRUE)
  }
  return(pwms)
}



#' Plot a genomic region surrounding a genomic variant, and potentially disrupted
#' motifs
#'
#' @param results The output of \code{motifbreakR}
#' @param rsid Character; the identifier (\code{SNP_id}) of the variant to be visualized
#' @param reverseMotif Logical; for motifs matched on the "-" strand, show the
#'   reverse complement of the motif (\code{TRUE}) or the motif reversed only
#'   (\code{FALSE})
#' @param effect Character; show motifs that are strongly affected \code{c("strong")},
#'   weakly affected \code{c("weak")}, or both \code{c("strong", "weak")}
#' @param altAllele Character; The default value of \code{NULL} uses the first (or only)
#'   alternative allele for the SNP to be plotted.
#' @seealso See \code{\link{motifbreakR}} for the function that produces output to be
#'   visualized here, also \code{\link{snps.from.rsid}} and \code{\link{snps.from.file}}
#'   for information about how to generate the input to \code{\link{motifbreakR}}
#'   function.
#' @details \code{plotMB} produces output showing the location of the variant on
#'   the chromosome, the reference and alternate sequence of the + strand, the
#'   footprint of any motif that is disrupted by the variant, and the DNA sequence
#'   motif(s), with the position of the variant marked on each motif.
#'   The \code{altAllele} argument is included for variants like rs1006140 where
#'   multiple alternate alleles exist, the reference allele is A, and the alternate
#'   can be G,T, or C. \code{plotMB} only plots one alternate allele at a time.
#' @return plots a figure representing the results of \code{motifbreakR} at the
#'   location of a single SNP, returns invisible \code{NULL}.
#' @examples
#' data(example.results)
#' example.results
#' \donttest{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' plotMB(results = example.results, rsid = "rs1006140", effect = "strong", altAllele = "C")
#' }
#' @importFrom motifStack DNAmotifAlignment colorset motifStack plotMotifLogo plotMotifLogoStack
#' @importClassesFrom motifStack pfm marker
#' @import grDevices
#' @importFrom grid gpar
#' @importFrom Gviz IdeogramTrack SequenceTrack GenomeAxisTrack HighlightTrack AnnotationTrack plotTracks
#' @export
plotMB <- function(results, rsid, reverseMotif = TRUE, effect = c("strong", "weak"), altAllele = NULL) {
  result <- results[results$SNP_id %in% rsid]
  if(is.null(altAllele)) {
    altAllele <- result$ALT[[1]]
  }
  res_filter <- (result$ALT == altAllele) & (result$effect %in% effect)
  result <- result[res_filter]

  motif.starts <- getAlleleSpecific(result, "motifCoordVariant", coord = "start", offset = -1L * start(result))
  motif.ends <- getAlleleSpecific(result, "motifCoordVariant", coord = "end", offset = -1L * (start(result) + getGaps(result)))

  result <- result[order(motif.starts, motif.ends)]

  chromosome <- as.character(seqnames(result))[[1]]
  genome.package <- attributes(result)$genome.package
  genome.bsgenome <- getExportedValue(genome.package, genome.package)

  seq.len <- max(length(result$ALT[[1]]), length(result$REF[[1]]))
  distance.to.edge <- 5 + seq.len
  from <- min(motif.starts) - distance.to.edge
  to <- max(motif.ends) + distance.to.edge + 1L
  pwmList <- attributes(result)$motifs
  pwm.names <- result$providerId
  results_motifs <- paste0(result$providerId, result$providerName)
  list_motifs <- paste0(mcols(pwmList)$providerId, mcols(pwmList)$providerName)
  pwms <- pwmList <- pwmList[match(results_motifs, list_motifs)]
  for (pwm.i in seq_along(pwms)) {
    pwm.name <- names(pwms[pwm.i])
    pwm.id <- mcols(pwms[pwm.name, ])$providerId
    pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
    doRev <- result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]
    doRevStrand <- ifelse(doRev$alleleDiff > 0,
                          doRev$motifStrandAlt,
                          doRev$motifStrandRef)
    doRev <- doRevStrand == 2L
    if (doRev) {
        pwm <- pwms[[pwm.i]]
        pwm <- pwm[, rev(1:ncol(pwm))]
        if (reverseMotif) {
          rownames(pwm) <- c("T", "G", "C", "A")
          pwm <- pwm[c("A", "C", "G", "T"), ]
          pwms[[pwm.i]] <- pwm
          names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], "-:rc")
        } else {
          pwms[[pwm.i]] <- pwm
          names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], "-:r")
        }
    }
  }
  pwms <- lapply(names(pwms), function(x, pwms=pwms) {new("pfm", mat = pwms[[x]],
                                                          name = x)}, pwms)
  pwms <- DNAmotifAlignment.2snp(pwms, result)
  pwmwide <- max(vapply(pwms, function(x) { ncol(x@mat)}, integer(1)))

  markerStart <- pintersect(getAlleleSpecific(result, "variantCoordMotif", coord = "seq"),
                            rep.int(IRanges(0, pwmwide), length(result)))
  markerEnd <- max(start(markerStart) + max(width(markerStart)) - 1L)
  markerStart <- max(start(markerStart))

  varType <- result$varType[[1]]
  varType <- switch(varType,
                    Deletion = "firebrick",
                    Insertion = "springgreen4",
                    SNV = "dodgerblue3",
                    Other = "gray13")
  markerRect <- new("marker", type = "rect",
                    start = markerStart + 1L,
                    stop = markerEnd + 1L,
                    gp = gpar(lty = 2,
                              fill = NA,
                              lwd = 3,
                              col = varType))
  for (pwm.i in seq_along(pwms)) {
    pwms[[pwm.i]]@markers <- list(markerRect)
  }
  g <- genome(genome.bsgenome)[[1]]

  ideoT <- try(IdeogramTrack(genome = g, chromosome = chromosome), silent = TRUE)
  if (inherits(ideoT, "try-error")) {
    backup.band <- data.frame(chrom = chromosome, chromStart = 0,
                              chromEnd = length(genome.bsgenome[[chromosome]]),
                              name = chromosome, gieStain = "gneg")
    ideoT <- IdeogramTrack(genome = g, chromosome = chromosome, bands = backup.band)
  }

  ### blank alt sequence
  altseq <- genome.bsgenome[[chromosome]]

  ### Replace longer sections
  at <- IRanges(start = start(result[1]), width = width(result[1]))

  axisT <- GenomeAxisTrack(exponent = 0)
  seqT <- SequenceTrack(genome.bsgenome,
                        fontcolor = colorset("DNA", "blindnessSafe"), add53 = TRUE,
                        chromosome = chromosome)
  sub_bases <- result$ALT[[1]]
  alt_diff <- 0L
  if(result$varType[[1]] == "Deletion") {
    alt_diff <- length(result$REF[[1]]) - length(result$ALT[[1]])
    addedN <- DNAString(paste0(rep.int(".", alt_diff), collapse = ""))
    sub_bases <- c(result$ALT[[1]], addedN)
    alt_diff <- 0L
  } else if (result$varType[[1]] == "Insertion") {
    alt_diff <- length(result$ALT[[1]]) - length(result$REF[[1]])
    addedN <- DNAString(paste0(rep.int("-", alt_diff), collapse = ""))
    refseq <- genome.bsgenome[[chromosome]]
    refseq <- DNAStringSet(replaceAt(x = refseq, at = at, c(result$REF[[1]], addedN)))
    names(refseq) <- chromosome
    rm(axisT)
    seqT <- SequenceTrack(refseq,
                          genome = g,
                          fontcolor = colorset("DNA", "blindnessSafe"),
                          chromosome = chromosome,
                          add53 = TRUE)
    names(seqT) <- "reference"

  }
  altseq <- DNAStringSet(replaceAt(x = altseq, at = at, sub_bases))

  names(altseq) <- chromosome
  seqAltT <- SequenceTrack(altseq, genome = g,
                           fontcolor = colorset("DNA", "blindnessSafe"),
                           chromosome = chromosome)
  hirange <- result[1]
  end(hirange) <- end(hirange) + alt_diff
  hiT <- HighlightTrack(trackList = list(seqT, seqAltT),
                        range = hirange)
  selectingfun <- selcor
  detailfun <- addPWM.stack

  motif_ids <- names(pwmList)
  names(motif_ids) <- mcols(pwmList)$providerName
  motif_ids <- motif_ids[result$providerName]
  presult <- result
  ranges(presult) <- IRanges(start = motif.starts, end = motif.ends)
  strand(presult) <- "*"
  pres_cols <- DataFrame(feature = paste(presult$geneSymbol, "motif", sep = "_"),
                         group = presult$providerName,
                         id = motif_ids)
  presult <- GRanges(seqnames = seqnames(presult[1]),
                     ranges = ranges(presult))
  mcols(presult) <- pres_cols

  motifT <- AnnotationTrack(presult,
                            fun = detailfun,
                            detailsFunArgs = list(pwm_stack = pwms),
                            name = result$SNP_id[[1]],
                            selectFun = selectingfun,
                            reverseStacking = FALSE,
                            stacking = "squish")

  if (exists("axisT")) {
    track_list <- list(ideoT, motifT, hiT, axisT)
  } else {
    track_list <- list(ideoT, motifT, hiT)
  }

  plotTracks(track_list, from = from, to = to, showBandId = TRUE,
             cex = 1, col.main = "darkgrey",
             labelpos = "below", chromosome = chromosome,
             fontcolor.item="black",
             collapse = FALSE, min.width = 1, featureAnnotation = "feature", cex.feature = 0.8,
             details.size = 0.85, detailsConnector.pch = NA, detailsConnector.lty = 0,
             shape = "box", cex.title = 1.1)
  return(invisible(NULL))
}

#' Run Shiny version of the motifbreakR package
#' @return returns a \code{\link[shiny]{shinyAppDir}} that launches the shiny app when printed.
#' @examples
#' library(motifbreakR)
#'
#' app <- shiny_motifbreakR()
#'
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
#' @import shiny bslib bsicons
#' @importFrom DT renderDT DTOutput datatable formatRound JS
#' @importFrom BSgenome available.genomes
#' @export
shiny_motifbreakR <- function() {
  appDir <- system.file("shiny", "shiny_motif", package = "motifbreakR")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `motifbreakR`.", call. = FALSE)
  }
  shinyAppDir(appDir)
}

#' Export motifbreakR results to csv or tsv
#'
#' @param results The output of \code{motifbreakR}
#' @param file Character; the file name of the destination file
#' @param format Character; one of tsv (tab separated values) or csv (comma separated values)
#' @return \code{exportMBtable} produces an output file containing the output
#' of the \code{motifbreakR} function.
#' @seealso See \code{\link{exportMBbed}} for the function that exports the
#' \code{motifbreakR} results as a BED file, colored by selected score.
#' @examples
#' data(example.results)
#' example.results
#' \donttest{
#' exportMBtable(example.results, file = tempfile(fileext = ".tsv"), format = "tsv")
#' }
#' @importFrom utils write.csv write.table
#' @export
exportMBtable <- function(results, file, format = "tsv") {
  if(missing(file)) {stop("select output file location")}
  format <- match.arg(format, c("csv", "tsv"))
  sep <- switch(format,
                csv = ",",
                tsv = "\t")
  results <- as.data.frame(results, row.names = NULL)
  results <- results[, !colnames(results) %in% "width"]
  results$start <- results$start - 1
  if("matchingCellType" %in% colnames(results)) {
    results$matchingBindingEvent <- vapply(results$matchingCellType, function(x) {
      collapsed_data <- lapply(x, function(y) {
        paste0(y, collapse = "; ")})
      collapsed_data <- paste0(names(collapsed_data), ":(", collapsed_data, ")")
      collapsed_data <- paste(collapsed_data, collapse = "; ")
      ifelse(collapsed_data == ":(NA)", NA_character_, collapsed_data)
    }, FUN.VALUE = character(1))
    results$matchingCellType <- NULL
  }
  if(format == "tsv") {
    write.table(x = results, file = file, quote = FALSE, sep = sep, row.names = FALSE)
  } else {
    write.csv(x = results, file = file, row.names = FALSE)
  }
}

get_color_values <- function(bed_score, color_set) {
  color_values <- quantile(bed_score, probs = seq(0,1,length.out = 9))
  cut_scores <- cut(bed_score, color_values, include.lowest = TRUE)
  cut_points <- levels(cut_scores)
  bed_colors <- factor(cut_scores,
                       labels = color_set)
  color_scale <- levels(bed_colors)
  names(color_scale) <- cut_points
  return(list(bed_colors = bed_colors, color_scale = color_scale))
}

#' Export motifbreakR variants to bed file
#'
#' @param results The output of \code{\link{motifbreakR}}
#' @param file Character; the file name of the destination file
#' @param name Character; name for the BED track, defaults to "motifbreakR results"
#' @param color Character; one of ref_sig (\code{pValueRef}), alt_sig
#' (\code{pValueAlt}), best_sig (lowest between \code{pValueRef} and
#' \code{pValueAlt}), (each of which require pre-computation of p-values with
#' \code{\link{calculatePvalue}}), or ref_score (\code{pctRef}), alt_score
#' (\code{pctAlt}), best_score (highest between \code{pctRef} and \code{pctAlt}),
#' or the default value of effect_size (\code{alleleDiff}).
#' @return \code{exportMBbed} produces an output BED file, with diverging color
#' scale for effect_size (blue representing stronger binding in \code{REF}, red
#' representing stronger binding in \code{ALT}), or a sequential color scale
#' otherwise (low values as purple, high values as yellow). The score column is
#' either the effect_size (\code{alleleDiff} column), the -log10(p-value)
#' (capped at 10), corresponding to \code{pValueRef}, \code{pValueAlt}, or the
#' best match of the two, or the score \code{pctRef}, \code{pctAlt}, or the
#' highest match of the two. The name column is formatted
#' \code{SNP_id:REF/ALT:providerId}. Additionally a color key is returned
#' indicating the range of values for each color output.
#' @seealso See \code{\link{exportMBtable}} for the function that exports the
#' full \code{motifbreakR} results as a tab or comma separated table file.
#' @examples
#' data(example.results)
#' example.results
#' \donttest{
#' exportMBbed(example.results, file = tempfile(fileext = ".bed"), color = "effect_size")
#' }
#' @export
exportMBbed <- function(results, file, name = NULL, color = "effect_size") {
  if(missing(file)) {stop("select output file location")}
  color <- match.arg(color, c("ref_sig", "alt_sig", "best_sig",
                              "ref_score", "alt_score", "best_score",
                              "effect_size"))
  sequential_colors <- c("#440154",
                         "#46337E",
                         "#365C8D",
                         "#277F8E",
                         "#1FA187",
                         "#4AC16D",
                         "#9FDA3A",
                         "#FDE725")
  diverging_colors <- c("#001260",
                        "#034A85",
                        "#4A8FB2",
                        "#BDD6E3",
                        "#E7C6B2",
                        "#C98157",
                        "#9D3709",
                        "#590007")
  bed_name <- paste(results$SNP_id, paste(results$REF, results$ALT, sep = "/"), sep = ":")
  bed_name <- paste(bed_name, results$providerId, sep = ":")
  results_bed <- results
  mcols(results_bed) <- NULL
  results_bed$name <- bed_name

  if(color == "effect_size") {
    bed_score <- results$alleleDiff
    color_values <- get_color_values(bed_score, diverging_colors)
  } else {
    if(color %in% c("ref_sig", "alt_sig", "best_sig")) {
      if(!("pValueRef" %in% names(mcols(results))))
        stop('incorrect results format; please rerun analysis with filterp=TRUE')
      if(any(is.na(results$pValueRef)))
        stop('run calculatePvalue before exporting data with p-values')
    }
    bed_score <- switch(color,
                        ref_sig = -log10(results$pValueRef),
                        alt_sig = -log10(results$pValueAlt),
                        best_sig = -log10(pmin(results$pValueRef, results$pValueAlt)),
                        ref_score = results$pctRef,
                        alt_score = results$pctAlt,
                        best_score = ifelse(results$pctAlt < results$pctRef,
                                            results$pctRef,
                                            results$pctAlt))
    color_values <- get_color_values(bed_score, sequential_colors)
  }
  results_bed$score <- bed_score
  results_bed$itemRgb <- color_values$bed_colors
  motif_sources <- levels(factor(results$dataSource))
  data_source <- unique(genome(results))
  results_trackline <- new("BasicTrackLine",
                           itemRgb = TRUE,
                           useScore = TRUE,
                           visibility = "pack",
                           name = ifelse(is.null(name), "motifbreakR results", name),
                           description = paste0("motifbreakR results colored by ", color, " in ",
                                                paste0(data_source, collapse = ", "),
                                                " from motif sources: ",
                                                paste0(motif_sources, collapse = ", ")))
  export.bed(results_bed, con = file, trackLine = results_trackline)
  return(color_values$color_scale)
}

#' Find Corresponding TF Binding From The ReMap2022 Project
#'
#' @param results The output of \code{motifbreakR}
#' @param genome Character; one of:
#' \code{hg38} or \code{hg19} for Homo sapiens,
#' \code{mm10} or \code{mm39} for Mus musculus,
#' \code{dm6} for Drosophila melanogaster,
#' \code{TAIR10_TF} or \code{TAIR10_HISTONE} for Arabidopsis thaliana
#' @param TFClass Logical;  The user may optionally query an expanded
#' motif/transcription factor relationship encompassing the entire potential
#' transcription factor family as implemented by \code{\link[MotifDb]{MotifDb}} based on
#' TFClass.
#' @details \code{TFClass} argument works for objects loaded in from the
#' \code{MotifDb} package. \code{hg19} and \code{mm39} are data from liftOver.
#'
#' The ReMap catalogues (2022, 2020, 2018, 2015) are under CC BY-NC 4.0
#' international license, as described in ReMap.
#'
#' The CC BY-NC 4.0 license correspond to the following terms:
#' Attribution — You must give appropriate credit, provide a link to the
#' license, and indicate if changes were made. You may do so in any reasonable
#' manner, but not in any way that suggests the licensor endorses you or your
#' use.
#' NonCommercial — You may not use the material for commercial purposes.
#' No additional restrictions — You may not apply legal terms or technological
#' measures that legally restrict others from doing anything the license
#' permits.
#'
#' @seealso \code{\link[MotifDb]{associateTranscriptionFactors}} for information about
#' TFClass. \url{https://remap.univ-amu.fr/} for details about ReMap2022.
#' @return the results GenomicRanges object output by \code{\link{motifbreakR}}
#' with the additional columns:
#'  \item{matchingBindingEvent}{The name of the transcription factor that binds
#'  over the motif, or \code{NA} if none}
#'  \item{matchingCellType}{A list corresponding in length to the number of
#'  transcription factors in \code{matchingBindingEvent} indicating the
#'  biotype/celltype that the transcription factor binding was found in.}
#' @examples
#' data(example.results)
#' \donttest{
#' example.results <- findSupportingRemapPeaks(example.results,
#'                                             genome = "hg19",
#'                                             TFClass = TRUE)
#' }
#' @export

findSupportingRemapPeaks <- function(results, genome, TFClass = FALSE) {
  genome <- match.arg(genome, c("hg38", "hg19", "mm10", "mm39", "dm6", "TAIR10_TF", "TAIR10_HISTONE"))
  switch(genome,
         hg19 = message("hg19 peaks have been lifted over from hg38"),
         mm39 = message("mm39 peaks have been lifted over from mm10"))

  remap_links <- list(
    hg38      = "https://remap.simoncoetzee.com/remap2022_nr_macs2_hg38_v1_0.bed.gz",
    hg19      = "https://remap.simoncoetzee.com/remap2022_nr_macs2_hg19_v1_0.bed.gz",
    mm10      = "https://remap.simoncoetzee.com/remap2022_nr_macs2_mm10_v1_0.bed.gz",
    mm39      = "https://remap.simoncoetzee.com/remap2022_nr_macs2_mm39_v1_0.bed.gz",
    dm6       = "https://remap.simoncoetzee.com/remap2022_nr_macs2_dm6_v1_0.bed.gz",
    TAIR10_TF = "https://remap.simoncoetzee.com/remap2022_nr_macs2_TAIR10_v1_0.bed.gz",
    TAIR10_HISTONE = "https://remap.simoncoetzee.com/remap2022_histone_nr_macs2_TAIR10_v1_0.bed.gz"
  )

  results$matchingBindingEvent <- NA
  results$matchingCellType <- NA

  uniqmb <- results
  mcols(uniqmb) <- NULL
  names(uniqmb) <- NULL
  uniqmb <- sort(unique(uniqmb))

  remap_binding <- loadPeakFile(remap_links, genome)
  remap_binding <- subsetByOverlaps(remap_binding, uniqmb)

  if(length(remap_binding) < 1) return(results)

  mb_overlaps <- findOverlaps(results, remap_binding)
  variant_to_peak <- unique(data.frame(variant = names(results[queryHits(mb_overlaps), ]), tf.binding = remap_binding[subjectHits(mb_overlaps), ]$name))
  tf_ct <- stringr::str_split(variant_to_peak$tf.binding, ":")
  variant_to_peak$tf.binding <- vapply(tf_ct, FUN = `[`, character(1), 1)
  variant_to_peak$tf.celltype <- vapply(tf_ct, FUN = `[`, character(1), 2)
  variant_to_peak$tf.celltype <- lapply(variant_to_peak$tf.celltype, function(x) {stringr::str_split(x, ",", simplify = T)[1, ]})
  mcol_variant_to_gene <- data.frame(variant = names(results), tf.gene = results$geneSymbol, motif = results$providerId)

  if(("manuallyCuratedGeneMotifAssociationTable" %in% slotNames(attributes(results)$motifs)) &
      TFClass &
      (nrow(attributes(results)$motifs@manuallyCuratedGeneMotifAssociationTable) > 1)) {
    mb_genelist <- attributes(results)$motifs@manuallyCuratedGeneMotifAssociationTable
    mb_genelist <- unique(mb_genelist[mb_genelist$motif %in% results$providerId, c("motif", "tf.gene")])
    mb_genelist_aug <- unique(data.frame(motif = results$providerId, tf.gene = results$geneSymbol))
    mb_genelist <- unique(rbind(mb_genelist, mb_genelist_aug))
    variant_to_gene <- data.frame(variant = names(results), motif = results$providerId)
    variant_to_gene_master <- merge(variant_to_gene, mb_genelist, all.x = TRUE, sort = FALSE)
    variant_to_gene <- unique(variant_to_gene_master[, c("variant", "tf.gene")])
  } else {
    variant_to_gene_master <- mcol_variant_to_gene
    variant_to_gene <- unique(variant_to_gene_master[, c("variant", "tf.gene")])
  }

  variant_mb_peak <- merge(variant_to_gene_master, variant_to_peak, all.x = TRUE, sort = FALSE)
  variant_mb_peak <- unique(variant_mb_peak[variant_mb_peak$tf.gene == variant_mb_peak$tf.binding, ])
  variant_mb_peak <- variant_mb_peak[!is.na(variant_mb_peak$variant), c("variant", "tf.binding", "tf.celltype", "motif")]

  split_vmbp <- split(variant_mb_peak, as.factor(variant_mb_peak$variant))

  for(variant in seq_along(split_vmbp)) {
    pvariant <- names(split_vmbp[variant])
    vgm <- split_vmbp[[variant]]
    mbv_sel <- which(mcol_variant_to_gene$variant == pvariant)
    vgm$motif <- factor(vgm$motif, levels = unique(mcol_variant_to_gene[mbv_sel, "motif"]))
    names(vgm$tf.celltype) <- vgm$tf.binding
    vgmm <- vgm[which(vgm$motif %in% mcol_variant_to_gene$motif), ]
    vgmm <- split(vgmm, vgmm$motif)
    vgmg <- vgm[which(vgm$tf.gene %in% mcol_variant_to_gene$tf.gene), ]
    vgmg <- split(vgmg, vgmg$motif)
    vgm <- mapply(rbind, vgmm, vgmg, SIMPLIFY = F)

    results[mbv_sel]$matchingBindingEvent <- lapply(vgm, function(x) {unique(x$tf.binding)})[results[mbv_sel]$providerId]
    results[mbv_sel]$matchingCellType <- lapply(vgm, function(x) {x$tf.celltype})[results[mbv_sel]$providerId]
    if(any(lengths(results[mbv_sel]$matchingBindingEvent) == 0)) {
      results[mbv_sel][lengths(results[mbv_sel]$matchingBindingEvent) == 0]$matchingBindingEvent <- NA
      results[mbv_sel][lengths(results[mbv_sel]$matchingCellType) == 0]$matchingCellType <- NA
    }
  }
  attributes(results)$peaks <- remap_binding
  return(results)
}

#' @importFrom tools R_user_dir
#' @importFrom BiocFileCache BiocFileCache
.get_cache <- function() {
  cache <- tools::R_user_dir("motifbreakR", which="cache")
  BiocFileCache(cache = cache, ask = F)
}

#' @importFrom BiocFileCache bfcquery bfcadd bfcneedsupdate bfcdownload bfcrpath bfcnew
cachePeakFile <- function(fileURL, genome) {
  bfc <- .get_cache()
  rname <- paste("remap2022", genome, sep = "_")
  rid <- bfcquery(bfc, rname, "rname")$rid
  if (!length(rid)) {
    rid <- names(bfcadd(bfc, rname, fileURL, ext = ".Rdata", download = FALSE))
  }
  needupdate <- tryCatch(bfcneedsupdate(bfc, rid),
                         error = function(cond) {
                           message("URL is down, using cache if available")
                           FALSE
                         })
  if (!isFALSE(needupdate)) {
    message("downloading peak file")
    bfcdownload(bfc, rid, ask = FALSE, FUN=convertPeakFile)
    message("peak download complete")
  }
  bfcrpath(bfc, rids = rid)
}

cacheMartObj <- function(genome) {
  bfc <- .get_cache()
  genome_ver <- ifelse(is.null(genome), "hg38", genome)
  rname <- paste("ensembl_mart", genome_ver, sep = "_")
  rid <- bfcquery(bfc, rname, "rname")$rid
  if (!length(rid)) {
    rid <- names(bfcnew(bfc, rname, ext = ".Rdata"))
    bmsnp <- useEnsembl(biomart = 'snps', dataset = "hsapiens_snp", version = genome)
    saveRDS(bmsnp, file = bfcrpath(bfc, rids = rid))
  }
  days.old <- as.numeric(Sys.Date() - as.Date(bfcquery(bfc, field = "rid", rid)$create_time))
  if (days.old > 30) {
    bmsnp <- useEnsembl(biomart = 'snps', dataset = "hsapiens_snp", version = genome)
    saveRDS(bmsnp, file = bfcrpath(bfc, rids = rid))
  }
  bfcrpath(bfc, rids = rid)
}

loadPeakFile <- function(url_list, genome) {
  peak_path <- cachePeakFile(url_list[[genome]], genome)
  remap_peaks <- readRDS(peak_path)
  return(remap_peaks)
}

loadMartObj <- function(genome) {
  mart_path <- cacheMartObj(genome)
  mart_obj <- readRDS(mart_path)
  return(mart_obj)
}

#' @importFrom vroom vroom
convertPeakFile <- function(from, to) {
  message("processing peak file")
  remap_peaks <- GRanges(vroom(from, delim = "\t",
                               col_names = c("chr", "start", "end", "name",
                                             "score", "strand", "tstart",
                                             "tend", "color"),
                               col_types = c("ciiciciic"),
                               col_select = c(chr, start, end, name),
                               progress = FALSE))
  saveRDS(remap_peaks, file = to)
  TRUE
}

