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

#' @importFrom matrixStats colRanges
#' @importFrom matrixStats colMaxs colMins
wScore <- function(snp.seq, ppm, offset, ppm.range = NULL, calcp=TRUE) {
  returnScore <- scoreMotif(snp.seq, ppm, ncol(ppm), offset = offset)
  if(calcp){
    return(sum(returnScore))
  } else {
    return((sum(returnScore)-ppm.range[1])/(ppm.range[2] - ppm.range[1]))
  }
}

scoreMotif.a <- function(snp.seq, ppm, len, offset = 1) {
  snp.seq <- snp.seq[offset:(offset + len - 1)]
  ## diag code
  position.probs <- c(ppm[snp.seq, ])[1 + 0L:(len - 1L) * (len + 1)]
  return(position.probs)
}
#' @importFrom compiler cmpfun
scoreMotif <- cmpfun(scoreMotif.a, options = list(optimize = 3))

## take a sequence and score all its windows for a pwm
scoreAllWindows <- function(snp.seq, snp.seq.rc, pwm,
                            from = "default", to = "default",
                            pwm.range = NULL, calcp=TRUE) {
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
    window.scores[i - from + 1] <- wScore(snp.seq, pwm, offset = i, ppm.range = pwm.range, calcp=calcp)
    window.scores.rc[i - from + 1] <- wScore(snp.seq.rc, pwm, offset = i, ppm.range = pwm.range, calcp=calcp)
  }
  all.window.scores <- matrix(data = c(window.scores, window.scores.rc), nrow = 2,
                              ncol = m, byrow = TRUE, dimnames = list(c("top", "bot"), from:to))
  return(all.window.scores)
}

## filter percent results to a threshold, e.g. 0.8 (%80), report only max value
## accept output of scoreWindows OR pctPwm

maxThresholdWindows <- function(window.frame, threshold = 0.8) {
  if (max(window.frame) > threshold) {
    arraymax <- which.max(window.frame)  ## the array index of max value
    arraycol <- floor(arraymax/2) + arraymax%%2  ## dereference the column (window coord)
    strand <- (2 * (arraymax%%2)) - 1  ## top or bottom strand from odd or even idx
    c(window = arraycol, strand = strand)
  } else {
    c(window = 0L, strand = 0L)
  }
}

## An evaluator function for SNP effect

snpEff <- function(allelR, allelA) {
  effect <- allelR - allelA
  if (abs(effect) < 0.4) {
    description <- "neut"
  } else {
    if (abs(effect) >= 0.7) {
      description <- "strong"
    } else {
      description <- "weak"
    }
  }
}

#' @import methods
#' @import GenomicRanges
#' @import S4Vectors
#' @import BiocGenerics
#' @import IRanges
#' @importFrom Biostrings getSeq replaceLetterAt reverseComplement complement
#' @importFrom TFMPvalue TFMpv2sc
scoreSnpList <- function(fsnplist, pwmList, method = "default", bkg = NULL,
                         threshold = 1e-3, show.neutral = FALSE, verbose = FALSE,
                         genome.bsgenome=NULL, pwmList.pc = NULL, pwmRanges = NULL, filterp=TRUE) {
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
    if(!filterp) {
      res.el$pctRef <- as.numeric(NA)
      res.el$pctAlt <- as.numeric(NA)
    } else {
      res.el$pctRef <- as.numeric(NA)
      res.el$pctAlt <- as.numeric(NA)
      res.el$scoreRef <- as.numeric(NA)
      res.el$scoreAlt <- as.numeric(NA)
      res.el$Refpvalue <- as.numeric(NA)
      res.el$Altpvalue <- as.numeric(NA)
    } 
    res.el$alleleRef <- as.numeric(NA)
    res.el$alleleAlt <- as.numeric(NA)
    res.el$effect <- as.character(NA)
    for (pwm.i in seq_along(pwmList)) {
      ## REF
      pwm.basic <- pwmList[[pwm.i]]
      pwm <- pwmList.pc[[pwm.i]]
      len <- ncol(pwm)
      thresh <- threshold[[pwm.i]]
      ref.windows <- scoreAllWindows(snp.ref, snp.ref.rc, pwm,
                                     from = k - len + 1, to = k, pwm.range = pwmRanges[[pwm.i]],
                                     calcp=filterp)
      hit.ref <- maxThresholdWindows(ref.windows, threshold = thresh)
      ## ALT
      alt.windows <- scoreAllWindows(snp.alt, snp.alt.rc, pwm,
                                     from = k - len + 1, to = k, pwm.range = pwmRanges[[pwm.i]],
                                     calcp=filterp)
      hit.alt <- maxThresholdWindows(alt.windows, threshold = thresh)
      if(hit.ref[["strand"]] != 0L && hit.alt[["strand"]] != 0L) {
        ref.strand <- (-hit.ref[["strand"]] + 3) / 2
        alt.strand <- (-hit.alt[["strand"]] + 3) / 2
        bigger <- ref.windows[ref.strand, hit.ref[["window"]]] > alt.windows[alt.strand, hit.alt[["window"]]]
        if(bigger) {
          hit <- hit.ref
        } else {
          hit <- hit.alt
        }
      } else {
        if (hit.ref[["strand"]] != 0L && hit.alt[["strand"]] == 0L) {
          hit <- hit.ref
        } else {
          if(hit.alt[["strand"]] != 0L && hit.ref[["strand"]] == 0L) {
            hit <- hit.alt
          } else {
            hit <-NULL
          }
        }
      }
      if(!is.null(hit)) {
        result <- res.el[pwm.i]
        snp.pos <- len - hit[["window"]] + 1
        uniquename <- paste(names(result), result$dataSource, result$providerName, result$providerId, sep = "%%")
        if (hit[["strand"]] == 1L) {
          allelR <- pwm.basic[as.character(result$REF), snp.pos]
          allelA <- pwm.basic[as.character(result$ALT), snp.pos]
        } else {
          allelR <- pwm.basic[as.character(complement(result$REF)), snp.pos]
          allelA <- pwm.basic[as.character(complement(result$ALT)), snp.pos]
        }
        effect <- snpEff(allelR, allelA)
        if (effect == "neut") {
          if (show.neutral) {
            res.el.e[[uniquename]] <- updateResults(result, snp.ref, snp.pos,
                                                    hit, ref.windows, alt.windows, allelR, allelA, effect, len,
                                                    k, pwm, calcp=filterp)
          }
        } else {
          res.el.e[[uniquename]] <- updateResults(result, snp.ref, snp.pos,
                                                hit, ref.windows, alt.windows, allelR, allelA, effect, len,
                                                k, pwm, calcp=filterp)
        }
      }
    }
  }
  resultSet <- unlist(GRangesList(as.list.environment(res.el.e)), use.names = FALSE)
  if (length(resultSet) < 1) {
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist),
                    "with 0 potentially disruptive matches to", length(unique(resultSet$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
    return(NULL)
  } else {
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
#' @importFrom TFMPvalue TFMsc2pv
updateResults <- function(result, snp.seq, snp.pos, hit, ref.windows, alt.windows,
                          allelR, allelA, effect, len, k, pwm, calcp) {
  strand.opt <- c("+", "-")
  strand <- (-hit[["strand"]] + 3) / 2
  strand(result) <- strand.opt[strand]
  mresult <- mcols(result)
  mresult[["snpPos"]] <- start(result)
  mresult[["motifPos"]] <- as.integer(snp.pos)
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
  if(calcp) {
    mresult[["scoreRef"]] <- ref.windows[strand, hit[["window"]]]
    mresult[["scoreAlt"]] <- alt.windows[strand, hit[["window"]]]
    mresult[["Refpvalue"]] <- NA
    mresult[["Altpvalue"]] <- NA
    pwmrange <- colSums(colRanges(pwm))
    mresult[["pctRef"]] <- (mresult[["scoreRef"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
    mresult[["pctAlt"]] <- (mresult[["scoreAlt"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
  } else {
    mresult[["pctRef"]] <- ref.windows[strand, hit[["window"]]]
    mresult[["pctAlt"]] <- alt.windows[strand, hit[["window"]]]
  }
  mresult[["alleleRef"]] <- allelR
  mresult[["alleleAlt"]] <- allelA
  mresult[["effect"]] <- effect
  mcols(result) <- mresult
  return(result)
}


#' Predict The Disruptiveness Of Single Nucleotide Polymorphisms On
#' Transcription Factor Binding Sites.
#'
#' @param snpList The output of \code{snps.from.rsid} or \code{snps.from.file}
#' @param pwmList An object of class \code{MotifList} containing the motifs that
#'   you wish to interrogate
#' @param threshold Numeric; the maximum p-value for a match to be called or a minimum score threshold
#' @param method Character; one of \code{default}, \code{log}, \code{ic}, or \code{notrans}; see
#'   details.
#' @param bkg Numeric Vector; the background probabilites of the nucleotides
#'   used with method=\code{log} method=\code{ic}
#' @param filterp Logical; filter by p-value instead of by pct score.
#' @param show.neutral Logical; include neutral changes in the output
#' @param verbose Logical; if running serially, show verbose messages
#' @param BPPARAM a BiocParallel object see \code{\link[BiocParallel]{register}}
#'   and see \code{getClass("BiocParallelParam")} for additional parameter
#'   classes.  Try \code{BiocParallel::registered()} to see what's availible and
#'   for example \code{BiocParallel::bpparam("SerialParam")} would allow serial
#'   evaluation.
#' @seealso See \code{\link{snps.from.rsid}} and \code{\link{snps.from.file}} for
#'   information about how to generate the input to this function and
#'   \code{\link{plotMB}} for information on how to visualize it's output
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
#' \eqn{s_{i} \in \{ A,T,C,G \}, i = 1,\ldots n}{s_i in {A,T,G,C}, i = 1 \ldots n}.
#' Each column of
#' \eqn{M} contains the frequencies of each letter in each position.
#'
#' Commonly in the literature sequences are scored as the sum of log
#' probabilities:
#'
#' \strong{Equation 1}
#'
#' \deqn{F( s,M ) = \sum_{i = 1}^{n}{\log( \frac{M_{s_{i},i}}{b_{s_{i}}} )}}{
#' F( s,M ) = \sum_(i = 1)^n log ((M_s_i,_i)/b_s_i)}
#'
#' where \eqn{b_{s_{i}}}{b_s_i} is the background frequency of letter \eqn{s_{i}}{s_i} in
#' the genome of interest. This method can be specified by the user as
#' \code{method='log'}.
#'
#' As an alternative to this method, we introduced a scoring method to
#' directly weight the score by the importance of the position within the
#' match sequence. This method of weighting is accessed by specifying
#' \code{method='ic'} (information content). A general representation
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
#' \deqn{p_{s} = ( M_{s_{i},i} ) \textrm{\ \ \ where\ \ \ } \frac{i = 1,\ldots n}{s_{i} \in \{ A,C,G,T \}}}{
#' p_s = ( M_s_i,_i ) where (i = 1 \ldots n)/(s_i in {A,C,G,T})}
#'
#' and second, for each \eqn{M} a constant vector of weights
#' \eqn{\omega_{M} = ( \omega_{1},\omega_{2},\ldots,\omega_{n} )}{\omega_M = ( \omega_1, \omega_2, \ldots, \omega_n)}.
#'
#' There are two methods for producing \eqn{\omega_{M}}{\omega_M}. The first, which we
#' call weighted sum, is the difference in the probabilities for the two
#' letters of the polymorphism (or variant), \emph{i.e.}
#' \eqn{\Delta p_{s_{i}}}{\Delta p_s_i}, or the difference of the maximum and minimum
#' values for each column of \eqn{M}:
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
#' \deqn{\omega_{i} = \sum_{j \in \{ A,C,G,T \}}^{}{M_{j,i}\log_2( \frac{M_{j,i}}{b_{i}} )}\textrm{\ \ \ \ \ where\ \ \ \ \ }i = 1,\ldots n}{
#' \omega_i = \sum_{j in {A,C,G,T}} {M_(j,i)} log2(M_(j,i)/b_i) where i = 1 \ldots n}
#'
#' where \eqn{b_{i}}{b_i} is again the background frequency of the letter \eqn{i}.
#'
#' Thus, there are 3 possible algorithms to apply via the \code{method}
#' argument. The first is the standard summation of log probabilities
#' (\code{method='log'}). The second and third are the weighted sum and
#' information content methods (\code{method='default'} and \code{method='ic'}) specified by
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
#' (\eqn{F( s_{\textsc{ref}},M )}{F( s_ref,M )} and
#' \eqn{F( s_{\textsc{alt}},M )}{F( s_alt,M )}), and provides the matrix scores
#' \eqn{p_{s_{\textsc{ref}}}}{p_s_ref} and \eqn{p_{s_{\textsc{alt}}}}{p_s_alt} of the SNP (or
#' variant). The scores are scaled as a fraction of scoring range 0-1 of
#' the motif matrix, \eqn{M}. If either of
#' \eqn{F( s_{\textsc{ref}},M )}{F( s_ref,M )} and
#' \eqn{F( s_{\textsc{alt}},M )}{F( s_alt,M )} is greater than a user-specified
#' threshold (default value of 0.85) the SNP is reported. By default
#' \pkg{motifbreakR} does not display neutral effects,
#' (\eqn{\Delta p_{i} < 0.4}{\Delta p_i < 0.4}) but this behaviour can be
#' overridden.
#'
#' Additionally, now, with the use of \code{\link{TFMPvalue-package}}, we may filter by p-value of the match.
#' This is unfortunately a two step process. First, by invoking \code{filterp=TRUE} and setting a threshold at
#' a desired p-value e.g 1e-4, we perform a rough filter on the results by rounding all values in the PWM to two
#' decimal place, and calculating a scoring threshold based upon that. The second step is to use the function \code{\link{calculatePvalue}()}
#' on a selection of results which will change the \code{Refpvalue} and \code{Altpvalue} columns in the output from \code{NA} to the p-value
#' calculated by \code{\link{TFMsc2pv}}.  This can be (although not always) a very memory and time intensive process if the algorithm doesn't converge rapidly.
#'
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
#'  \item{pctRef}{The score as determined by the scoring method, when the sequence contains the reference SNP allele, normalized to a scale from 0 - 1. If \code{filterp = FALSE},
#'  this is the value that is thresholded.}
#'  \item{pctAlt}{The score as determined by the scoring method, when the sequence contains the alternate SNP allele, normalized to a scale from 0 - 1. If \code{filterp = FALSE},
#'  this is the value that is thresholded.}
#'  \item{scoreRef}{The score as determined by the scoring method, when the sequence contains the reference SNP allele}
#'  \item{scoreAlt}{The score as determined by the scoring method, when the sequence contains the alternate SNP allele}
#'  \item{Refpvalue}{p-value for the match for the pctRef score, initially set to \code{NA}. see \code{\link{calculatePvalue}} for more information}
#'  \item{Altpvalue}{p-value for the match for the pctAlt score, initially set to \code{NA}. see \code{\link{calculatePvalue}} for more information}
#'  \item{alleleRef}{The proportional frequency of the reference allele at position \code{motifPos} in the motif}
#'  \item{alleleAlt}{The proportional frequency of the alternate allele at position \code{motifPos} in the motif}
#'  \item{effect}{one of weak, strong, or neutral indicating the strength of the effect.}
#'  each SNP in this object may be plotted with \code{\link{plotMB}}
#'  @examples
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
#' @importFrom BiocParallel bplapply
#' @importFrom stringr str_length str_trim
#' @export
motifbreakR <- function(snpList, pwmList, threshold=0.85, filterp = FALSE,
                        method = "default", show.neutral = FALSE, verbose = FALSE,
                        bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                        BPPARAM=bpparam()) {
  if(.Platform$OS.type == "windows" && inherits(BPPARAM, "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n",
            "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpworkers(BPPARAM)
  num.snps <- length(snpList)
  if(num.snps < cores) {
    cores <- num.snps
  }
  if(class(BPPARAM) == "SnowParam") {
    bpstart(BPPARAM)
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  genome.package <- attributes(snpList)$genome.package
  if(requireNamespace(eval(genome.package), quietly = TRUE, character.only = TRUE)) {
    genome.bsgenome <- eval(parse(text = paste(genome.package, genome.package, sep="::")))
  } else {
    stop(paste0(eval(genome.package), " is the genome selected for this snp list and \n",
                "is not present on your environment. Please load it and try again."))
  }
  snpList <- sapply(suppressWarnings(split(snpList, 1:cores)), list)
  bkg <- bkg[c('A', 'C', 'G', 'T')]

  scounts <- as.integer(mcols(pwmList)$sequenceCount)
  scounts[is.na(scounts)] <- 20L
  pwmList.pc <- Map(function(pwm, scount) {
                      pwm <- (pwm * scount + 0.25)/(scount + 1)
                    }, pwmList, scounts)
  if(method == "ic") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm, b=bkg) {
                          omegaic <- colSums(pwm * log2(pwm/b))
                        })
  }
  if(method == "default") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm) {
                          omegadefault <- colMaxs(pwm) - colMins(pwm)
                        })
  }
  if(method == "log") {
    pwmList.pc <- lapply(pwmList.pc, function(pwm, b) {
                           pwm <- log(pwm) - log(b)
                         }, b = bkg)
    pwmOmegas <- 1
  }
  if(method == "notrans") {
    pwmOmegas <- 1
  }
  pwmList.pc <- Map(function(pwm, omega) {
                      if(length(omega) == 1 && omega == 1) {
                        return(pwm)
                      } else {
                        omegamatrix <- matrix(rep(omega, 4), nrow = 4, byrow=TRUE)
                        pwm <- pwm * omegamatrix
                      }
                    }, pwmList.pc, pwmOmegas)
  if(filterp) {
    pwmRanges <- Map(function(pwm, omega) {
                       x <- colSums(colRanges(pwm))
                       return(x)
                     }, pwmList.pc, pwmOmegas)
    pwmList.pc2 <- lapply(pwmList.pc, round, digits = 2)
    pwmThresh <- lapply(pwmList.pc2, TFMpv2sc, pvalue = threshold, bg = bkg, type = "PWM")
    pwmThresh <- Map("+", pwmThresh, -0.02)
  } else {
    pwmRanges <- Map(function(pwm, omega) {
                       x <- colSums(colRanges(pwm))
                       return(x)
                     }, pwmList.pc, pwmOmegas)
    pwmThresh <- rep.int(threshold, times = length(pwmRanges))
  }
  pwmList@listData <- lapply(pwmList, function(pwm) { pwm <- rbind(pwm, N=0); return(pwm) })
  pwmList.pc <- lapply(pwmList.pc, function(pwm) { pwm <- rbind(pwm, N=0); return(pwm) })
  x <- try(bplapply(snpList, scoreSnpList, pwmList = pwmList, threshold = pwmThresh,
                    method = method, bkg = bkg, show.neutral = show.neutral, pwmList.pc = pwmList.pc,
                    verbose = ifelse(cores == 1, verbose, FALSE), genome.bsgenome = genome.bsgenome,
                    pwmRanges = pwmRanges, filterp = filterp, BPPARAM=BPPARAM))
  if(inherits(x, "try-error")) {
    if(class(BPPARAM) == "SnowParam") {
      bpstop(BPPARAM)
    }
    stop(attributes(x)$condition)
  }
  if(class(BPPARAM) == "SnowParam") {
	  bpstop(BPPARAM)
  }
  drops <- sapply(x, is.null)
  x <- x[!drops]
  pwmList@listData <- lapply(pwmList, function(pwm) { pwm <- pwm[c("A", "C", "G", "T"), ]; return(pwm) })
  pwmList.pc <- lapply(pwmList.pc, function(pwm) { pwm <- pwm[c("A", "C", "G", "T"), ]; return(pwm) })
  if (length(x) > 1) {
    x <- unlist(GRangesList(unname(x)))
    x <- x[order(match(names(x), names(snpList)), x$geneSymbol), ]
    attributes(x)$genome.package <- genome.package
    attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% unique(x$providerId) &
                                      mcols(pwmList)$providerName %in% unique(x$providerName), ]
    attributes(x)$scoremotifs <- pwmList.pc[names(attributes(x)$motifs)]
  } else {
    if (length(x) == 1L) {
      x <- x[[1]]
      attributes(x)$genome.package <- genome.package
      attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% unique(x$providerId) &
                                        mcols(pwmList)$providerName %in% unique(x$providerName), ]
      attributes(x)$scoremotifs <- pwmList.pc[names(attributes(x)$motifs)]
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


#' Calculate the significance of the matches for the reference and alternate alleles for the for their PWM
#'
#' @param results The output of \code{motifbreakR} that was run with \code{filterp=TRUE}
#' @param background Numeric Vector; the background probabilites of the nucleotides
#' @return a GRanges object. The same Granges object that was input as \code{results}, but with
#'  \code{Refpvalue} and \code{Altpvalue} columns in the output modified from \code{NA} to the p-value
#'  calculated by \code{\link{TFMsc2pv}}.
#' @seealso See \code{\link{TFMsc2pv}} from the \pkg{TFMPvalue} package for
#'   information about how the p-values are calculated.
#' @details This function is intended to be used on a selection of results produced by \code{\link{motifbreakR}}, and
#' this can be (although not always) a very memory and time intensive process if the algorithm doesn't converge rapidly.
#' @source H{\'e}l{\`e}ne Touzet and Jean-St{\'e}phane Varr{\'e} (2007) Efficient and accurate P-value computation for Position Weight Matrices.
#'  Algorithms for Molecular Biology, \bold{2: 15}.
#' @examples
#' data(example.results)
#' rs2661839 <- example.results[names(example.results) %in% "rs2661839"]
#' rs2661839 <- calculatePvalue(rs2661839)
#'
#' @export
calculatePvalue <- function(results,
                            background = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)) {
  if(!("scoreRef" %in% names(mcols(results)))) {
    stop('incorrect results format; please rerun analysis with filterp=TRUE')
  } else {
    pwmListmeta <- mcols(attributes(results)$motifs, use.names=TRUE)
    pwmList <- attributes(results)$scoremotifs
    pvalues <- lapply(results, function(result, pwmList, pwmListmeta, bkg) {
                        pwm.id <- result$providerId
                        pwm.name.f <- result$providerName
                        pwmmeta <- pwmListmeta[pwmListmeta$providerId == pwm.id & pwmListmeta$providerName == pwm.name.f, ]
                        pwm <- pwmList[[rownames(pwmmeta)[1]]]
                        ref <- TFMsc2pv(pwm, mcols(result)[["scoreRef"]], bg = bkg, type="PWM")
                        alt <- TFMsc2pv(pwm, mcols(result)[["scoreAlt"]], bg = bkg, type="PWM")
                        return(data.frame(ref=ref, alt=alt))
    }, pwmList=pwmList, pwmListmeta=pwmListmeta, bkg = background)
    pvalues.df <- base::do.call("rbind", c(pvalues, make.row.names = FALSE))
    results$Refpvalue <- pvalues.df[, "ref"]
    results$Altpvalue <- pvalues.df[, "alt"]
    return(results)
  }
}

#addPWM <- function(identifier, ...) {
#  motif.name <- identifier
#  ps.file <- paste(tempdir(), motif.name, sep = "/")
#  PostScriptTrace(ps.file, paste0(ps.file, ".xml"))
#  motif.figure <- readPicture(paste0(ps.file, ".xml"))
#  grid.picture(motif.figure[-1], distort = FALSE)
#}


#' @importFrom grImport PostScriptTrace readPicture grid.picture
addPWM.stack <- function(identifier, index, GdObject, ...) {
  psloc <- tempdir()
  ps.file <- paste(psloc, "stack.ps", sep = "/")
  PostScriptTrace(ps.file, paste0(ps.file, ".xml"))
  motif.figure <- readPicture(paste0(ps.file, ".xml"))
  motif.figure <- motif.figure[-1]
  snppos <- sapply(sapply(mcols(GdObject@range)[, "feature"], strsplit, "@"), "[[", 2)
  motif.i <- 1
  highlight <- snppos[[motif.i]]
  paths <- motif.figure@paths
  for(path.i in seq_along(paths)) {
    if (inherits(paths[[path.i]], "PictureText")) {
      if (paths[[path.i]]@string == highlight) {
        path <- paths[[path.i]]
        path@rgb <- "#FF0000"
        path@letters <- lapply(path@letters, function(letters) {
          letters@rgb <- "#FF0000"
          return(letters)
        })
        paths[[path.i]] <- path
        motif.i <- motif.i + 1
        if(motif.i <= length(snppos)) {
          highlight <- snppos[[motif.i]]
        } else {
          motif.i <- length(snppos)
        }
      }
    }
  }
  motif.figure@paths <- paths
  grid.picture(motif.figure, distort = FALSE)
}

#  grid.rect(gp = gpar(fill = rgb(1,0,0,1)), draw = T, width = 70, height = 1000, default.units = "native", x = 2000, y=0)
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
#  pos <- mcols(result)$motifPos
#  pos <- pos[as.logical(strand(result) == "+")][1]
#  browser()
  for(pwm.i in seq_along(pwms)) {
    pwm <- pwms[[pwm.i]]@mat
    ## get pwm info from result data
    pwm.name <- pwms[[pwm.i]]@name
    pwm.name <- str_replace(pwm.name, pattern = "-:rc$", replacement = "")
    pwm.name <- str_replace(pwm.name, pattern = "-:r$", replacement = "")
    pwm.info <- attributes(result)$motifs
    pwm.id <- mcols(pwm.info[pwm.name, ])$providerId
    pwm.name <- mcols(pwm.info[pwm.name, ])$providerName
    start.offset <- start(result[result$providerId == pwm.id & result$providerName == pwm.name, ]) - from
    end.offset <- to - end(result[result$providerId == pwm.id & result$providerName == pwm.name, ])
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
#' @param effect Character; show motifs that are strongly effected \code{c("strong")},
#'   weakly effected \code{c("weak")}, or both \code{c("strong", "weak")}
#' @seealso See \code{\link{motifbreakR}} for the function that produces output to be
#'   visualized here, also \code{\link{snps.from.rsid}} and \code{\link{snps.from.file}}
#'   for information about how to generate the input to \code{\link{motifbreakR}}
#'   function.
#' @details \code{plotMB} produces output showing the location of the SNP on the
#'   chromosome, the surrounding sequence of the + strand, the footprint of any
#'   motif that is disrupted by the SNP or SNV, and the DNA sequence motif(s)
#' @return plots a figure representing the results of \code{motifbreakR} at the
#'   location of a single SNP, returns invisible \code{NULL}.
#' @examples
#' data(example.results)
#' example.results
#' \dontrun{
#' plotMB(example.results, "rs2661839", effect = "strong")
#' }
#' @import motifStack
#' @import grDevices
#' @importFrom Gviz IdeogramTrack SequenceTrack GenomeAxisTrack HighlightTrack
#'   AnnotationTrack plotTracks
#' @export
plotMB <- function(results, rsid, reverseMotif = TRUE, effect = c("strong", "weak")) {
  g <- genome(results)[[1]]
  result <- results[names(results) %in% rsid]
  stackmotif <- TRUE
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
    getmotifs <- mcols(pwmList)$providerId %in% result$providerId & mcols(pwmList)$providerName %in% result$providerName
    pwms <- pwmList[getmotifs, ]
    pwms <- pwms[order(match(paste0(mcols(pwms)$providerId, mcols(pwms)$providerName), paste0(result$providerId, result$providerName)))]
    if(reverseMotif) {
      for(pwm.i in seq_along(pwms)) {
        pwm.name <- names(pwms[pwm.i])
        pwm.id <- mcols(pwms[pwm.name, ])$providerId
        pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
        doRev <- as.logical(strand(result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]) == "-")
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
        pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
        doRev <- as.logical(strand(result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]) == "-")
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
    pwmwide <- max(sapply(pwms, function(x) { ncol(x@mat)}))
    psloc <- tempdir()
    postscript(paste(psloc, "stack.ps", sep = "/"), width = pwmwide * (1/3), height = 2 * length(pwm.names), paper="special", horizontal = FALSE)
    plotMotifLogoStack.2(pwms, ncex=1.0)
    dev.off()
  } else {
    getmotifs <- mcols(pwmList)$providerId %in% result$providerId & mcols(pwmList)$providerName %in% result$providerName
    pwmList <- pwmList[getmotifs, ]
    pwmList <- pwmList[order(match(mcols(pwmList)$providerId, result$providerId))]
    for(pwm.i in seq_along(pwmList)) {
      pwm <- pwmList[pwm.i]
      if (reverseMotif) {
        pwm.id <- mcols(pwm)$providerId
        pwm.name.f <- mcols(pwm)$providerName
        doRev <- as.logical(strand(result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]) == "-")
        if(doRev) {
          pwm.mat <- pwm[[1]]
          pwm.mat <- pwm.mat[, rev(1:ncol(pwm.mat))]
          rownames(pwm.mat) <- c("T", "G", "C", "A")
          pwm.mat <- pwm.mat[c("A", "C", "G", "T"), ]
          pwm[[pwm.i]] <- pwm.mat
          names(pwm)[1] <- paste0(names(pwm)[1], "-:rc")
          pwm@elementMetadata@rownames <- names(pwm)[1]
          pwmList@elementMetadata[pwm.i,] <- pwm@elementMetadata
          pwmList@listData[pwm.i] <- pwm@listData
          names(pwmList@listData)[pwm.i] <- names(pwm)[1]
        }
      } else {
        pwm.id <- mcols(pwm)$providerId
        pwm.name.f <- mcols(pwm)$providerName
        doRev <- as.logical(strand(result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]) == "-")
        if(doRev) {
          pwm.mat <- pwm[[1]]
          pwm.mat <- pwm.mat[, rev(1:ncol(pwm.mat))]
          pwm[[pwm.i]] <- pwm.mat
          names(pwm)[1] <- paste0(names(pwm)[1], "-:r")
          pwm@elementMetadata@rownames <- names(pwm)[1]
          pwmList@elementMetadata[pwm.i,] <- pwm@elementMetadata
          pwmList@listData[pwm.i] <- pwm@listData
          names(pwmList@listData)[pwm.i] <- names(pwm)[1]
        }
      }
      p <- new("pfm", mat = pwm[[1]], name = names(pwm[1]))
      #message(paste(temp.dir, names(pwm[1]), sep="/"))
      psloc <- tempdir()
      pwmwide <- ncol(p@mat)
      postscript(paste(psloc, "stack.ps", sep = "/"), width = pwmwide * (7/11), height =  4, paper="special", horizontal = FALSE)
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
  ### blank alt sequence
  altseq <- genome.bsgenome[[chromosome]]
  wherereplace <- rep.int(TRUE, length(altseq))
  wherereplace[result$snpPos[[1]]] <- FALSE
  altseq <- replaceLetterAt(altseq, at = wherereplace, letter = rep.int("N", sum(wherereplace)))
  altseq <- replaceLetterAt(altseq, at = !wherereplace, letter = result$ALT[[1]])
  altseq <- DNAStringSet(altseq)
  names(altseq) <- chromosome
  seqAltT <- SequenceTrack(altseq, fontcolor = c(colorset("DNA", "auto"), N="#FFFFFF"), chromosome = chromosome)
  axisT <- GenomeAxisTrack(exponent = 0)
  hiT <- HighlightTrack(trackList = list(seqT, seqAltT), start = result$snpPos[[1]], end = result$snpPos[[1]],
                        chromosome = chromosome)
  if(stackmotif){
    selectingfun <- selcor
    detailfun <- addPWM.stack
  } else {
    detailfun <- addPWM.stack
    selectingfun <- selall
  }
  getmotifs <- mcols(pwmList)$providerId %in% result$providerId & mcols(pwmList)$providerName %in% result$providerName
  motifT <- AnnotationTrack(result, id = names(pwmList)[getmotifs],
                            fun = detailfun, group = result$providerName,
                            feature = paste0("snp@",
                                             ifelse(strand(result) == "-",
                                                    str_length(str_trim(result$seqMatch)) - result$motifPos + 1,
                                                    result$motifPos)),
                            name = names(result)[[1]], selectFun = selectingfun)
  plotTracks(list(ideoT, motifT, hiT, axisT), from = from, to = to, showBandId = TRUE,
             cex.main = 0.8, col.main = "darkgrey",
             add53 = TRUE, labelpos = "below", chromosome = chromosome, groupAnnotation = "group",
             collapse = FALSE, min.width = 1, featureAnnotation = "feature", cex.feature = 0.8,
             details.size = ifelse(stackmotif, 0.85, 0.5), detailsConnector.pch = NA, detailsConnector.lty = ifelse(stackmotif, 0, 3),
             shape = "box", cex.group = 0.8, fonts = c("sans", "Helvetica"))
  return(invisible(NULL))
}

