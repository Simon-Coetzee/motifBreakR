## scoreMotif.R this module will score an arbitrary list of SNPs for motif
## disruption using a buffet of menu options for some sets of transcription
## factors

## optional: produce random SNP sets of equivalent size to the index SNP list to
## simulate background for enrichment calculations.

## define utility functions 'defaultOmega' (and related fxns) and maxPwm and
## minPwm

## reverse complement any string letter DNA
revcom <- function(ltr) {
  rclist <- list(A = "T", C = "G", G = "C", T = "A")
  rclist[[ltr]]
}

## omega weighting function creates weights to be passed into scoreMotif with pwm
## should be calculated only once per funcisnp run, per motif
defaultOmega <- function(ppm) {
  l <- ncol(ppm)
  omega <- rep(0, length(l))
  for (i in 1:ncol(ppm)) {
    omega[i] <- max(ppm[, i]) - min(ppm[, i])
  }
  omega
}

## information content, see Stormo 2000 PMID: 10812473 eq 2
defaultIC <- function(ppm, bkg) {
  ## for ppm with i columns containing letters a,c,g,t:
  IC <- colSums(ppm * log2(ppm/bkg + 1e-07))
  return(IC)
}

## maxPwm will return the max function of pwm given omega constraints like the
## defaultOmega fxn, it should be run at a high level w/in the script to reduce
## overhead
limitPwm <- function(pwm, limitFun, method = "default", bkg = NULL) {
  if (is.null(bkg)) {
    bkg <- c(0.25, 0.25, 0.25, 0.25)
  } else {
    bkg <- bkg[c(1, 2, 2, 1)]
  }
  names(bkg) <- c("A", "C", "G", "T")
  if (method != "default") {
    pwm <- pwm/bkg
  }
  PWMlimit <- apply(pwm, 2, limitFun)
  if (method != "default") {
    PWMlimit <- PWMlimit + 1e-07
  }
  if (method == "default") {
    omega <- defaultOmega(pwm)
    vscores <- PWMlimit * omega
  } else if (method == "log") {
    vscores <- log(PWMlimit)
  } else if (method == "IC") {
    vscores <- PWMlimit * defaultIC(pwm, bkg)
  }
  sum(vscores)
}


pctScore <- function(score, min.score, max.score) {
  (score - min.score)/(max.score - min.score)
}

## motif scorer, using MotifDb objects returns vector of probabilities based on
## the PPM
#' @importFrom compiler cmpfun
scoreMotif.a <- function(snp.seq, ppm, len, method = "default",
                       ag = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                       offset = 1) {
  snp.seq <- snp.seq[offset:(offset + len - 1)]
  ## diag code
  position.probs <- c(ppm[snp.seq, ])[1 + 0L:(len - 1L) * (len + 1)]
  return(position.probs)
}
scoreMotif <- cmpfun(scoreMotif.a, options = list(optimize = 3))

## create an all purpose scoring algorithm 'w' for weighted
wScore <- function(snp.seq, ppm, pwm.omega, offset, method = "default", bkg = NULL,
                   len) {
  if (is.null(bkg)) {
    bkg <- c(0.25, 0.25, 0.25, 0.25)
  } else {
    bkg <- bkg[c(1, 2, 2, 1)]
  }
  if(method != "default") {
    ppm <- ppm/bkg
  }
  vscores <- scoreMotif(snp.seq, ppm, len, method = method, ag = bkg, offset = offset)
  if (method != "default") {
    vscores <- vscores + 1e-07
  }
  if (method == "default") {
    returnScore <- sum(vscores * pwm.omega)
  }
  if (method == "log") {
    returnScore <- sum(log(vscores))
  }
  if (method == "IC") {
    returnScore <- sum(vscores * defaultIC(ppm, bkg))
  }
  return(returnScore)
}

## take a sequence and score all its windows for a pwm

scoreAllWindows <- function(snp.seq, snp.seq.rc, pwm, pwm.omega, pwm.bot, pwm.min,
                            method = "default", bkg = NULL, from = "default", to = "default") {
  ## frequently used variables;
  l <- ncol(pwm)
  if (from == "default") {
    from <- 1
  }
  ## if ( to=='default') { to <- max(sapply(position.matches, max)) - l }
  if (to == "default") {
    to <- from + l - 1
  }
  m <- to - from + 1  ## number of windows of width l
  ## define a temporary pair of vectors to store scores
  window.scores <- rep(NA, m)
  window.scores.rc <- rep(NA, m)
  for (i in from:to) {
    window.scores[i - from + 1] <- wScore(snp.seq, pwm, pwm.omega, offset = i,
                                          method, bkg, len = l)
    window.scores.rc[i - from + 1] <- wScore(snp.seq.rc, pwm, pwm.omega, offset = i,
                                             method, bkg, len = l)
  }
  all.window.scores <- matrix(data = c(window.scores, window.scores.rc), nrow = 2,
                              ncol = m, byrow = T, dimnames = list(c("top", "bot"), from:to))
  all.window.scores <- (all.window.scores - pwm.min)/pwm.bot
  return(all.window.scores)
}

## argument window.frame is the matrix return value from scoreAllWindows convert
## scores to percentage score

pctPwm <- function(window.frame, pwm.max, pwm.min, method = "default", bkg = NULL) {
  window.pcts <- pctScore(window.frame, pwm.min, pwm.max)
  window.pcts
}

## filter percent results to a threshold, e.g. 0.8 (%80), report only max value
## accept output of scoreWindows OR pctPwm

maxThresholdWindows <- function(window.frame, threshold = 0.8) {
  if (max(window.frame) > threshold) {
    arraymax <- which.max(window.frame)  ## the array index of max value
    arraycol <- floor(arraymax/2) + arraymax%%2  ## dereference the column (window coord)
    strand <- (2 * (arraymax%%2)) - 1  ## top or bottom strand from odd or even idx
    c(window = arraycol, strand = strand)
  } else c(window = 0L, strand = 0L)
}

## An evaluator function for SNP effect

snpEff <- function(allelR, allelA) {
  effect <- allelR - allelA
  if (abs(effect) < 0.4) {
    description <- "neut"
  } else {
    if (abs(effect) >= 0.7) {
      description <- "strg"
    } else {
      description <- "weak"
    }
  }
}

#' @import methods
#' @importFrom Biostrings getSeq replaceLetterAt reverseComplement
#' @importFrom GenomicRanges promoters GRangesList
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocGenerics width
#' @importFrom IRanges reverse
scoreSnpList <- function(fsnplist, pwmList, method = "default", bkg = NULL,
                         threshold = 0.9, show.neutral = FALSE, verbose = FALSE, genome.bsgenome=NULL) {
  if (!(Reduce("&", width(fsnplist$REF) == 1)) || !(Reduce("&", width(fsnplist$ALT) == 1))) {
    drop.snp.i <- !(width(fsnplist$REF) == 1) & (width(fsnplist$ALT) == 1)
    drop.snp <- names(fsnplist[!drop.snp.i])
    fsnplist <- fsnplist[drop.snp.i]
    warning(paste(drop.snp, "not valid for motif analysis...skipping", sep = " "))
  }
  k <- max(sapply(pwmList, ncol))
  ## Get Reference Sequence with k bp flanking
  snp.sequence.ref <- getSeq(genome.bsgenome, promoters(fsnplist, upstream = k - 1,
                                                        downstream = k))
  ## check that reference matches ref genome
  equals.ref <- getSeq(genome.bsgenome, fsnplist) == fsnplist$REF
  if (!Reduce("&", equals.ref)) {
    stop(paste(names(fsnplist[!equals.ref]), "reference allele does not match value in reference genome",
               sep = " "))
  }
  at <- matrix(FALSE, nrow = length(snp.sequence.ref), ncol = (k * 2) - 1)
  at[, k] <- TRUE
  snp.sequence.alt <- replaceLetterAt(snp.sequence.ref, at, fsnplist$ALT)
  snp.sequence.ref.rc <- lapply(snp.sequence.ref, reverseComplement)
  snp.sequence.alt.rc <- lapply(snp.sequence.alt, reverseComplement)
  snp.sequence.ref <- lapply(snp.sequence.ref, function(x) {
    strsplit(as.character(x), "")[[1]]
  })
  snp.sequence.ref.rc <- lapply(snp.sequence.ref.rc, function(x) {
    strsplit(as.character(x), "")[[1]]
  })
  snp.sequence.alt <- lapply(snp.sequence.alt, function(x) {
    strsplit(as.character(x), "")[[1]]
  })
  snp.sequence.alt.rc <- lapply(snp.sequence.alt.rc, function(x) {
    strsplit(as.character(x), "")[[1]]
  })
  resultSet <- NULL
  resultSet <- lapply(fsnplist, function(x, l = length(pwmList)) {
    return(rep(x, l))
  })
  pwmList.omega <- lapply(pwmList, defaultOmega)
  pwmList.max <- lapply(pwmList, limitPwm, limitFun = max, method = method, bkg = bkg)
  pwmList.min <- lapply(pwmList, limitPwm, limitFun = min, method = method, bkg = bkg)
  pwmList.bot <- mapply("-", pwmList.max, pwmList.min)
  strand.opt <- c("+", "-")
  #set.sig <- rep(NA, length(resultSet))
  res.el.e <- new.env()
  for (snp.map.i in seq_along(snp.sequence.alt)) {
    snp.ref <- snp.sequence.ref[[snp.map.i]]
    snp.alt <- snp.sequence.alt[[snp.map.i]]
    snp.ref.rc <- snp.sequence.ref.rc[[snp.map.i]]
    snp.alt.rc <- snp.sequence.alt.rc[[snp.map.i]]
    res.el <- resultSet[[snp.map.i]]
    res.el$snpPos <- as.integer(NA)
    res.el$motifPos <- as.integer(NA)
    res.el$motifID <- mcols(pwmList)$providerID
    res.el$geneSymbol <- mcols(pwmList)$geneSymbol
    res.el$dataSource <- mcols(pwmList)$dataSource
    res.el$providerName <- mcols(pwmList)$providerName
    res.el$providerId <- mcols(pwmList)$providerId
    res.el$seqMatch <- as.character(NA)
    res.el$pctRef <- as.numeric(NA)
    res.el$pctAlt <- as.numeric(NA)
    res.el$alleleRef <- as.numeric(NA)
    res.el$alleleAlt <- as.numeric(NA)
    res.el$effect <- as.character(NA)
    #res.el.e <- list2env(as.list(split(res.el, 1:length(res.el))), parent = emptyenv())
    #pwm.sig <- rep(NA, length(res.el))
    for (pwm.i in seq_along(pwmList)) {
      ## REF
      pwm <- pwmList[[pwm.i]]
      len <- ncol(pwm)
      ref.windows <- scoreAllWindows(snp.ref, snp.ref.rc, pwm, pwmList.omega[[pwm.i]],
                                     pwmList.bot[[pwm.i]], pwmList.min[[pwm.i]], method, bkg = bkg, from = k -
                                       len + 1, to = k)
      hit.ref <- maxThresholdWindows(ref.windows, threshold = threshold)
      ## ALT
      alt.windows <- scoreAllWindows(snp.alt, snp.alt.rc, pwm, pwmList.omega[[pwm.i]],
                                     pwmList.bot[[pwm.i]], pwmList.min[[pwm.i]], method, bkg = bkg, from = k -
                                       len + 1, to = k)
      hit.alt <- maxThresholdWindows(alt.windows, threshold = threshold)
      if (hit.ref[["strand"]] == 0L && hit.alt[["strand"]] == 0L) {
       # pwm.sig[[pwm.i]] <- FALSE
      } else {
        if (hit.ref[["strand"]] != 0L) {
          hit <- hit.ref
          snp.pos <- len - hit[["window"]] + 1
          result <- res.el[pwm.i]
          uniquename <- paste(names(result), result$dataSource, result$providerName, result$providerId, sep = "%%")
          if (hit[["strand"]] == 1L) {
            allelR <- pwm[as.character(result$REF), snp.pos]
            allelA <- pwm[as.character(result$ALT), snp.pos]
          } else {
            allelR <- pwm[as.character(complement(result$REF)), snp.pos]
            allelA <- pwm[as.character(complement(result$ALT)), snp.pos]
          }
          effect <- snpEff(allelR, allelA)
          if (effect == "neut") {
            if (show.neutral) {
              res.el.e[[uniquename]] <- updateResults(result, snp.ref, snp.pos,
                                                      hit, ref.windows, alt.windows, allelR, allelA, effect, len,
                                                      k, pwm)
              #res.el.tmp <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
              #                            hit, ref.windows, alt.windows, allelR, allelA, effect, len,
              #                            k, pwm)
              #mcols(res.el)[pwm.i, ] <- mcols(res.el.tmp)
              #ranges(res.el)[pwm.i] <- ranges(res.el.tmp)
              #strand(res.el)[pwm.i] <- strand(res.el.tmp)
              #res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
              #                               hit, ref.windows, alt.windows, allelR, allelA, effect, len,
              #                               k, pwm)
             # pwm.sig[[pwm.i]] <- TRUE
            } else {
             # pwm.sig[[pwm.i]] <- FALSE
            }
          } else {
            res.el.e[[uniquename]] <- updateResults(result, snp.ref, snp.pos,
                                                    hit, ref.windows, alt.windows, allelR, allelA, effect, len,
                                                    k, pwm)
            #             res.el.tmp <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
            #                                         hit, ref.windows, alt.windows, allelR, allelA, effect, len,
            #                                         k, pwm)
            #             mcols(res.el)[pwm.i, ] <- mcols(res.el.tmp)
            #             ranges(res.el)[pwm.i] <- ranges(res.el.tmp)
            #             strand(res.el)[pwm.i] <- strand(res.el.tmp)
            #res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
            #                               hit, ref.windows, alt.windows, allelR, allelA, effect, len,
            #                               k, pwm)
            #pwm.sig[[pwm.i]] <- TRUE
          }
        }
        if (hit.alt[["strand"]] != 0L && hit.alt[["window"]] != hit.ref[["window"]]) {
          hit <- hit.alt
          snp.pos <- len - hit[["window"]] + 1
          result <- res.el[pwm.i]
          uniquename <- paste(names(result), result$dataSource, result$providerName, result$providerId, sep = "%%")
          if (hit[["strand"]] == 1L) {
            allelR <- pwm[as.character(result$REF), snp.pos]
            allelA <- pwm[as.character(result$ALT), snp.pos]
          } else {
            allelR <- pwm[as.character(complement(result$REF)), snp.pos]
            allelA <- pwm[as.character(complement(result$ALT)), snp.pos]
          }
          effect <- snpEff(allelR, allelA)
          if (effect == "neut") {
            if (show.neutral) {
              res.el.e[[uniquename]] <- updateResults(result, snp.ref, snp.pos,
                                                      hit, ref.windows, alt.windows, allelR, allelA, effect, len,
                                                      k, pwm)
              #               res.el.tmp <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
              #                                           hit, ref.windows, alt.windows, allelR, allelA, effect, len,
              #                                           k, pwm)
              #               mcols(res.el)[pwm.i, ] <- mcols(res.el.tmp)
              #               ranges(res.el)[pwm.i] <- ranges(res.el.tmp)
              #               strand(res.el)[pwm.i] <- strand(res.el.tmp)
              #res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
              #                               hit, ref.windows, alt.windows, allelR, allelA, effect, len,
              #                               k, pwm)
              #pwm.sig[[pwm.i]] <- TRUE
            } else {
              #pwm.sig[[pwm.i]] <- FALSE
            }
          } else {
            res.el.e[[uniquename]] <- updateResults(result, snp.ref, snp.pos,
                                                    hit, ref.windows, alt.windows, allelR, allelA, effect, len,
                                                    k, pwm)
            #             res.el.tmp <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
            #                                         hit, ref.windows, alt.windows, allelR, allelA, effect, len,
            #                                         k, pwm)
            #             mcols(res.el)[pwm.i, ] <- mcols(res.el.tmp)
            #             ranges(res.el)[pwm.i] <- ranges(res.el.tmp)
            #             strand(res.el)[pwm.i] <- strand(res.el.tmp)
            #res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos,
            #                               hit, ref.windows, alt.windows, allelR, allelA, effect, len,
            #                               k, pwm)
            #pwm.sig[[pwm.i]] <- TRUE
          }
        }
      }
    }
#    if (verbose) {
#      message(names(res.el[1]))
#      message(paste0("number of broken motifs = ", sum(pwm.sig)))
#    }
    #res.el.sldfkjds <<- res.el.e
    #stop('here')

    #res.el <- res.el[pwm.sig]
    #if (length(res.el) > 0) {
    #  resultSet[[snp.map.i]] <- res.el
    #  set.sig[[snp.map.i]] <- TRUE
    #} else {
    #  set.sig[[snp.map.i]] <- FALSE
    #}
  }
  #resultSet <- resultSet[set.sig]
  resultSet <- unlist(GRangesList(as.list.environment(res.el.e)), use.names = F)
  if (length(resultSet) < 1) {
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist),
                    "with 0 potentially disruptive matches to", length(unique(resultSet$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
    return(NULL)
  } else {
    #resultSet <- unlist(GRangesList(resultSet), use.names = F)
    mcols(resultSet) <- mcols(resultSet)[, -(1:2)]
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist),
                    "with", length(resultSet), "potentially disruptive matches to", length(unique(resultSet$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
    return(resultSet)
  }
}

#' @importFrom stringr str_pad
#' @importFrom BiocGenerics start end start<- end<-
updateResults <- function(result, snp.seq, snp.pos, hit, ref.windows, alt.windows,
                          allelR, allelA, effect, len, k, pwm) {
  strand.opt <- c("+", "-")
  strand <- (-hit[["strand"]] + 3) / 2
  strand(result) <- strand.opt[strand]
  mresult <- mcols(result)
  mresult[["snpPos"]] <- start(result)
  mresult[["motifPos"]] <- snp.pos
  if (strand == 2) {
    matchs <- snp.seq[(k - hit[["window"]] + 1):(k - hit[["window"]] + len)]
    matchs[-hit[["window"]]] <- tolower(matchs[-hit[["window"]]])
    matchs <- paste(matchs, collapse = "")
    mresult[["seqMatch"]] <- str_pad(matchs, width = k + hit[["window"]] - 1, side = "right")
    ## pwm pos instead of snp pos
    start(result) <- start(result) - hit[["window"]] + 1
    end(result) <- end(result) - hit[["window"]] + len
  } else {
    matchs <- snp.seq[(k - len + hit[["window"]]):(k + hit[["window"]] - 1)]
    matchs[-snp.pos] <- tolower(matchs[-snp.pos])
    matchs <- paste(matchs, collapse = "")
    mresult[["seqMatch"]] <- str_pad(matchs, width = k + snp.pos - 1, side = "right")
    start(result) <- start(result) - snp.pos + 1
    end(result) <- end(result) - snp.pos + len
  }
  mresult[["pctRef"]] <- ref.windows[strand, hit[["window"]]]
  mresult[["pctAlt"]] <- alt.windows[strand, hit[["window"]]]
  mresult[["alleleRef"]] <- allelR
  mresult[["alleleAlt"]] <- allelA
  mresult[["effect"]] <- effect
  mcols(result) <- mresult
  return(result)
}


#' Predict The Disruptiveness Of Single Nucleotide Polymorphisms On
#' Transcription Factor Binding Sites.
#'
#' @param snpList The output of \code{snps.from.rsid} or \code{snps.from.bed}
#' @param pwmList An object of class \code{MotifList} containing the motifs that
#'   you wish to interrogate
#' @param threshold Numeric; the minimum disruptiveness score for which to
#'   report results
#' @param method Character; one of \code{default}, \code{log}, or \code{IC}; see
#'   details.
#' @param bkg Numeric Vector; the background probabilites of the nucleotides
#'   used with method=\code{log} method=\code{IC}
#' @param show.neutral Logical; include neutral changes in the output
#' @param verbose Logical; if running serially, show verbose messages
#' @param BPPARAM a BiocParallel object see \code{\link[BiocParallel]{register}}
#'   and see \code{getClass("BiocParallelParam")} for additional parameter
#'   classes.  Try \code{BiocParallel::registered()} to see what's availible and
#'   for example \code{BiocParallel::bpparam("SerialParam")} would allow serial
#'   evaluation.
#' @seealso See \code{\link{snps.from.rsid}} and \code{\link{snps.from.bed}} for
#'   information about how to generate the input to this function and \code{\link{MBplot}}
#'   for information on how to visualize it's output
#' @details \pkg{motifbreakR} works with position probability matrices (PPM). PPM
#' are derived as the fractional occurrence of nucleotides A,C,G, and T at
#' each position of a position frequency matrix (PFM). PFM are simply the
#' tally of each nucleotide at each position across a set of aligned
#' sequences. With a PPM, one can generate probabilities based on the
#' genome, or more practically, create any number of position specific
#' scoring matrices (PSSM) based on the principle that the PPM contains
#' information about the likelihood of observing a particular nucleotide at
#' a particular position of a true transcription factor binding site. What
#' follows is a discussion of the three different algorithms that may be
#' employed in calls to the \pkg{motifbreakR} function via the \code{method}
#' argument.
#'
#' Suppose we have a frequency matrix \eqn{M} of width \eqn{n} (\emph{i.e.} a
#' PPM as described above). Furthermore, we have a sequence \eqn{s} also of
#' length \eqn{n}, such that
#' \eqn{s_{i} \in \{ A,T,C,G \}, i = 1,\ldots n}. Each column of
#' \eqn{M} contains the frequencies of each letter in each position.
#'
#' Commonly in the literature sequences are scored as the sum of log
#' probabilities:
#'
#' \strong{Equation 1}
#'
#' \deqn{F( s,M ) = \sum_{i = 1}^{n}{\log( \frac{M_{s_{i},i}}{b_{s_{i}}} )}}
#'
#' where \eqn{b_{s_{i}}} is the background frequency of letter \eqn{s_{i}} in
#' the genome of interest. This method can be specified by the user as
#' \code{method='log'}.
#'
#' As an alternative to this method, we introduced a scoring method to
#' directly weight the score by the importance of the position within the
#' match sequence. This method of weighting is accessed by specifying
#' \code{method='IC'} (information content). A general representation
#' of this scoring method is given by:
#'
#' \strong{Equation 2}
#'
#' \deqn{F( s,M ) = p( s ) \cdot \omega_{M}}
#'
#' where \eqn{p_{s}} is the scoring vector derived from sequence \eqn{s} and matrix
#' \eqn{M}, and \eqn{w_{M}} is a weight vector derived from \eqn{M}. First, we
#' compute the scoring vector of position scores \eqn{p}:
#'
#' \strong{Equation 3}
#'
#' \deqn{p( s ) = ( M_{s_{i},i} ) \textrm{\ \ \ where\ \ \ } \frac{i = 1,\ldots n}{s_{i} \in \{ A,C,G,T \}}}
#'
#' and second, for each \eqn{M} a constant vector of weights
#' \eqn{\omega_{M} = ( \omega_{1},\omega_{2},\ldots,\omega_{n} )}.
#'
#' There are two methods for producing \eqn{\omega_{M}}. The first, which we
#' call weighted sum, is the difference in the probabilities for the two
#' letters of the polymorphism (or variant), \emph{i.e.}
#' \eqn{\Delta p_{s_{i}}}, or the difference of the maximum and minimum
#' values for each column of \eqn{M}:
#'
#' \strong{Equation 4.1}
#'
#' \deqn{\omega_{i} = \max \{ M_{i} \} - \min \{ M_{i} \}\textrm{\ \ \ \ where\ \ \ \ \ \ }i = 1,\ldots n}
#'
#' The second variation of this theme is to weight by relative entropy.
#' Thus the relative entropy weight for each column \eqn{i} of the matrix is
#' given by:
#'
#' \strong{Equation 4.2}
#'
#' \deqn{\omega_{i} = \sum_{j \in \{ A,C,G,T \}}^{}{M_{j,i}\log( \frac{M_{j,i}}{b_{i}} )}\textrm{\ \ \ \ \ where\ \ \ \ \ }i = 1,\ldots n}
#'
#' where \eqn{b_{i}} is again the background frequency of the letter \eqn{i}.
#'
#' Thus, there are 3 possible algorithms to apply via the \code{method}
#' argument. The first is the standard summation of log probabilities
#' (\code{method='log'}). The second and third are the weighted sum and
#' information content methods (\code{method=default} and \code{method='IC'}) specified by
#' equations 4.1 and 4.2, respectively. \pkg{motifbreakR} assumes a
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
#' (\eqn{F( s_{\textsc{ref}},M )} and
#' \eqn{F( s_{\textsc{alt}},M )}), and provides the matrix scores
#' \eqn{p_{s_{\textsc{ref}}}} and \eqn{p_{s_{\textsc{alt}}}} of the SNP (or
#' variant). The scores are scaled as a fraction of scoring range 0-1 of
#' the motif matrix, \eqn{M}. If either of
#' \eqn{F( s_{\textsc{ref}},M )} and
#' \eqn{F( s_{\textsc{alt}},M )} is greater than a user-specified
#' threshold (default value of 0.85) the SNP is reported. By default
#' \pkg{motifbreakR} displays only strong effects
#' (\eqn{\Delta p_{i} > 0.7}) but this behaviour can be
#' overridden.
#' @return a GRanges object containing:
#'  \item{REF}{the reference allele for the SNP}
#'  \item{ALT}{the alternate allele for the SNP}
#'  \item{snpPos}{the coordinates of the SNP}
#'  \item{motifPos}{the coordinates of the SNP within the TF binding motif}
#'  \item{geneSymbol}{the geneSymbol corresponding to the TF of the TF binding motif}
#'  \item{dataSource}{the source of the TF binding motif}
#'  \item{providerName, providerId}{the name and id provided by the source}
#'  \item{seqMatch}{the sequence on the 5' -> 3' direction of the "+" strand
#'  that corresponds to DNA at the position that the TF binding motif was found.}
#'  \item{pctRef}{The score as determined by the scoring method, when the sequence contains the reference SNP allele}
#'  \item{pctAlt}{The score as determined by the scoring method, when the sequence contains the alternate SNP allele}
#'  \item{alleleRef}{The proportional frequency of the reference allele at position \code{motifPos} in the motif}
#'  \item{alleleAlt}{The proportional frequency of the alternate allele at position \code{motifPos} in the motif}
#'  \item{effect}{one of weak, strong, or neutral indicating the strength of the effect.}
#'  each SNP in this object may be plotted with \code{\link{MBplot}}
#'  @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  library(SNPlocs.Hsapiens.dbSNP.20120608)
#'  # prepare variants
#'  data(pca.enhancer.snps) # loads snps.mb
#'  pca.enhancer.snps <- sample(snps.mb, 20)
#'  # Get motifs to interrogate
#'  data(hocomoco)
#'  motifs <- sample(hocomoco, 50)
#'  # run motifbreakR
#'  results <- motifbreakR(pca.enhancer.snps, motifs, threshold = 0.9, method = "IC")
#' @import BiocParallel
#' @importFrom BiocParallel bplapply
#' @importFrom BiocGenerics unlist sapply
#' @importFrom stringr str_length str_trim
#' @export
motifbreakR <- function(snpList, pwmList, threshold, method = "default",
                        bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                        show.neutral = FALSE, verbose = FALSE, BPPARAM=bpparam()) {
  if(.Platform$OS.type == "windows" && inherits(BPPARAM, "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n",
            "Windows, please supply an alternative BPPARAM"))
  }
  if(inherits(BPPARAM, "SnowParam")) {
  }
  cores <- bpworkers(BPPARAM)
  num.snps <- length(snpList)
  if(num.snps < cores) {
    cores <- num.snps
  }
  if(inherits(BPPARAM, "SnowParam")) {
    bpstart(BPPARAM)
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  genome.package <- attributes(snpList)$genome.package
  if(requireNamespace(eval(genome.package), quietly = TRUE, character.only = TRUE)) {
    genome.bsgenome <- eval(parse(text = paste(genome.package, genome.package, sep="::")))
  } else {
    stop(paste0(eval(genome.package), " is the genome selected for this snp list and \n",
                "is not present on your system. Please install and try again."))
  }
  snpList <- sapply(suppressWarnings(split(snpList, 1:cores)), list)
  bkg <- bkg[c('A', 'C', 'G', 'T')]
  x <- bplapply(snpList, scoreSnpList, pwmList = pwmList, threshold = threshold,
                method = method, bkg = bkg, show.neutral = show.neutral,
                verbose = ifelse(cores == 1, verbose, FALSE), genome.bsgenome = genome.bsgenome,
                BPPARAM=BPPARAM)
  if(inherits(BPPARAM, "SnowParam")) {
    bpstop(BPPARAM)
  }
  x.glb <<- x
  drops <- sapply(x, is.null)
  x <- x[!drops]
  if (length(x) > 1) {
    x <- unlist(GRangesList(unname(x)))
    x <- x[order(match(names(x), names(snpList)), x$geneSymbol), ]
    attributes(x)$genome.package <- genome.package
    attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% unique(x$providerId), ]
  } else {
    if (length(x) == 1L) {
      x <- x[[1]]
      attributes(x)$genome.package <- genome.package
      attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% unique(x$providerId), ]
    } else {
      warning("No SNP/Motif Interactions reached threshold")
      x <- NULL
    }
  }
  if (verbose && cores > 1) {
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

#' @importFrom grImport PostScriptTrace readPicture grid.picture
addPWM <- function(identifier, ...) {
  motif.name <- identifier
  ps.file <- paste(tempdir(), motif.name, sep = "/")
  PostScriptTrace(ps.file, paste0(ps.file, ".xml"))
  motif.figure <- readPicture(paste0(ps.file, ".xml"))
  grid.picture(motif.figure[-1], distort = FALSE)
}

addPWM.stack <- function(identifier, ...) {
  ps.file <- "/Users/coetzeesg/stack.ps"
  PostScriptTrace(ps.file, paste0(ps.file, ".xml"))
  motif.figure <- readPicture(paste0(ps.file, ".xml"))
  grid.picture(motif.figure[-1], distort = FALSE)
}

selcor <- function(identifier, GdObject, ... ) {
  if(identifier == mcols(GdObject@range)$id[[1]]) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


selall <- function(identifier, GdObject, ... ) {
    return(TRUE)
}

plotMotifLogoStack.2 <- function(pfms, ...) {
  pfms <- rev(pfms)
  n <- length(pfms)
  lapply(pfms, function(.ele) {
    if (class(.ele) != "pfm")
      stop("pfms must be a list of class pfm")
  })
  opar <- par(mfrow = c(n, 1), mar = c(3.5, 3.5, 1.5, 0.5))
  assign("tmp_motifStack_symbolsCache", list(), pos = ".GlobalEnv")
  for (i in 1:(n - 1)) {
    motifStack::plotMotifLogo(pfms[[n - i + 1]], motifName = pfms[[n - i + 1]]@name,
                              p=rep(0.25, 4), xlab = NA)
  }
  motifStack::plotMotifLogo(pfms[[1]], motifName = pfms[[1]]@name, p=rep(0.25, 4))
  rm(list = "tmp_motifStack_symbolsCache", pos = ".GlobalEnv")
  par(opar)
}

#' @importFrom stringr str_replace
DNAmotifAlignment.2snp <- function(pwms, result) {
  from <- min(start(result))
  to <- max(end(result))
  for(pwm.i in seq_along(pwms)) {
    pwm <- pwms[[pwm.i]]@mat
    ## get pwm info from result data
    pwm.name <- pwms[[pwm.i]]@name
    pwm.name <- str_replace(pwm.name, pattern = "-:rc$", replace = "")
    pwm.name <- str_replace(pwm.name, pattern = "-:r$", replace = "")
    pwm.info <- attributes(result)$motifs
    pwm.id <- mcols(pwm.info[pwm.name, ])$providerId
    start.offset <- start(result[result$providerId == pwm.id, ]) - from
    end.offset <- to - end(result[result$providerId == pwm.id, ])
    if(start.offset > 0){
      pwm <- cbind(matrix(c(0.25, 0.25, 0.25, 0.25), ncol = start.offset, nrow = 4), pwm)
    }
    if(end.offset > 0) {
      pwm <- cbind(pwm, matrix(c(0.25, 0.25, 0.25, 0.25), ncol = end.offset, nrow = 4))
    }
    pwms[[pwm.i]]@mat <- pwm
  }
  return(pwms)
}



#' Plot a genomic region surrounding a genomic variant, and potentially disrupted
#' motifs
#'
#' @param results The output of \code{motifbreakR}
#' @param rsid Character; the identifier of the variant to be visualized
#' @param reverseMotif Logical; if the motif is on the "-" strand show the
#'   the motifs as reversed \code{FALSE} or reverse complement \code{TRUE}
#' @param stackmotif Logical; show motifs vertically, aligned on the variant \code{TRUE}
#'   or horizontally.
#' @param effect Character; show motifs that are strongly effected \code{c("strg")},
#'   weakly effected \code{c("weak")}, or both \code{c("strg", "weak")}
#' @seealso See \code{\link{motifbreakR}} for the function that produces output to be
#'   visualized here, also \code{\link{snps.from.rsid}} and \code{\link{snps.from.bed}}
#'   for information about how to generate the input to \code{\link{motifbreakR}}
#'   function.
#' @details \code{plotMB} produces output showing the location of the SNP on the
#'   chromosome, the surrounding sequence of the + strand, the footprint of any
#'   motif that is disrupted by the SNP or SNV, and the DNA sequence motif(s)
#' @return plots a figure representing the results of \code{motifbreakR} at the
#'   location of a single SNP, returns \code{NULL}.
#' @importFrom grDevices postscript
#' @import motifStack
#' @importFrom Gviz IdeogramTrack SequenceTrack GenomeAxisTrack HighlightTrack
#'   AnnotationTrack plotTracks
#' @importFrom BiocGenerics strand
#' @export
plotMB <- function(results, rsid, reverseMotif = TRUE, stackmotif = FALSE, effect = c("strg", "weak")) {
  g <- genome(results)[[1]]
  result <- results[names(results) %in% rsid]
  if(length(result) <= 1) stackmotif <- FALSE
  result <- result[order(start(result), end(result)), ]
  result <- result[result$effect %in% effect]
  chromosome <- as.character(seqnames(result))[[1]]
  genome.package <- attributes(result)$genome.package
  genome.bsgenome <- eval(parse(text = genome.package))
  distance.to.edge <- max(result[1, ]$snpPos - min(start(result)), max(end(result)) - result[1, ]$snpPos) + 4
  from <- result[1, ]$snpPos - distance.to.edge + 1
  to <- result[1, ]$snpPos + distance.to.edge
  temp.dir <- tempdir()
  pwmList <- attributes(result)$motifs
  pwm.names <- result$providerId
  if(stackmotif && length(result) > 1) {
    pwms <- pwmList[mcols(pwmList)$providerId %in% pwm.names, ]
    pwms <- pwms[order(match(mcols(pwms)$providerId, result$providerId))]
    if(reverseMotif) {
      for(pwm.i in seq_along(pwms)) {
        pwm.name <- names(pwms[pwm.i])
        pwm.id <- mcols(pwms[pwm.name, ])$providerId
        doRev <- as.logical(strand(result[result$providerId == pwm.id, ]) == "-")
        if(doRev) {
          pwm <- pwms[[pwm.i]]
          pwm <- pwm[, rev(1:ncol(pwm))]
          rownames(pwm) <- c("T", "G", "C", "A")
          pwm <- pwm[c("A", "C", "G", "T"), ]
          pwms[[pwm.i]] <- pwm
          names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], "-:rc")
        }
      }
    } else {
      for(pwm.i in seq_along(pwms)) {
        pwm.name <- names(pwms[pwm.i])
        pwm.id <- mcols(pwms[pwm.name, ])$providerId
        doRev <- as.logical(strand(result[result$providerId == pwm.id, ]) == "-")
        if(doRev) {
          pwm <- pwms[[pwm.i]]
          pwm <- pwm[, rev(1:ncol(pwm))]
          pwms[[pwm.i]] <- pwm
          names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], "-:r")
        }
      }
    }
    pwms <- lapply(names(pwms), function(x, pwms=pwms) {new("pfm", mat=pwms[[x]],
                                                          name=x)}, pwms)
    pwms <- DNAmotifAlignment.2snp(pwms, result)
    postscript(paste("~", "stack.ps", sep = "/"), width = 7, height =  2*length(pwm.names), paper="special", horizontal = FALSE)
    plotMotifLogoStack.2(pwms, ncex=1.0)
    dev.off()
  } else {
    pwmList <- pwmList[mcols(pwmList)$providerId %in% pwm.names, ]
    pwmList <- pwmList[order(match(mcols(pwmList)$providerId, result$providerId))]
    for(pwm.i in seq_along(pwmList)) {
      pwm <- pwmList[pwm.i]
      if (reverseMotif) {
        pwm.id <- mcols(pwm)$providerId
        doRev <- as.logical(strand(result[result$providerId == pwm.id, ]) == "-")
        if(doRev) {
          pwm.mat <- pwm[[1]]
          pwm.mat <- pwm.mat[, rev(1:ncol(pwm.mat))]
          rownames(pwm.mat) <- c("T", "G", "C", "A")
          pwm.mat <- pwm.mat[c("A", "C", "G", "T"), ]
          pwmList[[pwm.i]] <- pwm.mat
          names(pwm)[1] <- paste0(names(pwm)[1], "-:rc")
          names(pwmList)[pwm.i] <- names(pwm)[1]
        }
      } else {
        pwm.id <- mcols(pwm)$providerId
        doRev <- as.logical(strand(result[result$providerId == pwm.id, ]) == "-")
        if(doRev) {
          pwm.mat <- pwm[[1]]
          pwm.mat <- pwm.mat[, rev(1:ncol(pwm.mat))]
          pwmList[[pwm.i]] <- pwm.mat
          names(pwm)[1] <- paste0(names(pwm)[1], "-:r")
          names(pwmList)[pwm.i] <- names(pwm)[1]
        }
      }
      p <- new("pfm", mat = pwm[[1]], name = names(pwm[1]))
      message(paste(temp.dir, names(pwm[1]), sep="/"))
      postscript(paste(temp.dir, names(pwm[1]), sep = "/"), width = 10, height = 3, horizontal = FALSE,
                 fonts = c("sans", "Helvetica"))
      motifStack::plotMotifLogo(p, motifName = p@name)
      dev.off()
    }
  }
  ideoT <- try(IdeogramTrack(genome = g, chromosome = chromosome), silent = TRUE)
  if(inherits(ideoT, "try-error")) {
    backup.band <- data.frame(chrom = chromosome, chromStart = 0, chromEnd = length(genome.bsgenome[[chromosome]]), name = chromosome, gieStain = "gneg")
    ideoT <- IdeogramTrack(genome = g, chromosome = chromosome, bands = backup.band)
  }
  seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", "auto"))
  axisT <- GenomeAxisTrack(exponent = 0)
  hiT <- HighlightTrack(seqT, start = result$snpPos[[1]], end = result$snpPos[[1]],
                        chromosome = chromosome)
  if(stackmotif){
    selectingfun <- selcor
    detailfun <- addPWM.stack
  } else {
    detailfun <- addPWM
    selectingfun <- selall
  }
  motifT <- AnnotationTrack(result, id = names(pwmList)[mcols(pwmList)$providerId %in% result$providerId],
                            fun = detailfun, group = result$providerId,
                            feature = paste0("snp@",
                                             ifelse(strand(result) == "-",
                                                    str_length(str_trim(result$seqMatch)) - result$motifPos + 1,
                                                    result$motifPos)),
                            name = paste0(names(result)[[1]], " ", result$REF[[1]], "/", result$ALT[[1]]), selectFun = selectingfun)
  plotTracks(list(ideoT, motifT, hiT, axisT), from = from, to = to, showBandId = TRUE,
             add53 = TRUE, labelpos = "below", chromosome = chromosome, groupAnnotation = "group",
             collapse = FALSE, min.width = 1, featureAnnotation = "feature", cex.feature = 0.8,
             details.size = ifelse(stackmotif, 0.9, 0.5), detailsConnector.pch = NA, detailsConnector.lty = ifelse(stackmotif, 0, 3),
             shape = "box", cex.group = 0.8, fonts = c("sans", "Helvetica"))
  return(invisible(NULL))
}



