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
   position.probs <- c(ppm[snp.seq, ])[1L + 0L:(len - 1L) * (len + 1L)]
   return(position.probs)
 }
 #' @importFrom compiler cmpfun
 scoreMotif <- cmpfun(scoreMotif.a, options = list(optimize = 3))

# take a sequence and score all its windows for a pwm
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

# filter percent results to a threshold, e.g. 0.8 (%80), report only max value
# accept output of scoreWindows OR pctPwm

 # maxThresholdWindows <- function(window.frame, threshold = 0.8) {
 #   if (max(window.frame) > threshold) {
 #     arraymax <- which.max(window.frame)  ## the array index of max value
 #     arraycol <- floor(arraymax/2) + arraymax%%2  ## dereference the column (window coord)
 #     strand <- (2 * (arraymax%%2)) - 1  ## top or bottom strand from odd or even idx
 #     c(window = arraycol, strand = strand)
 #   } else {
 #     c(window = 0L, strand = 0L)
 #   }
 # }

## An evaluator function for SNP effect

snpEff <- function(allelR, allelA) {
  score <- allelR - allelA
  if (abs(score) < 0.4) {
    return(list(score = score, effect = "neut"))
  } else {
    if (abs(score) >= 0.7) {
      return(list(score = score, effect = "strong"))
    } else {
      return(list(score = score, effect = "weak"))
    }
  }
}

scoreIndel <- function(pwm,
                       ref.seq, alt.seq,
                       hit.ref, hit.alt) {
  ref.windows <- scoreSeqWindows(ppm = pwm, seq = ref.seq)
  alt.windows <- scoreSeqWindows(ppm = pwm, seq = alt.seq)
  score <- alt.windows[hit.alt$strand, hit.alt$window] - ref.windows[hit.ref$strand, hit.ref$window]
  if (abs(score) >= 0.7) {
    return(list(score = score, effect = "strong"))
  } else if (abs(score) > 0.4) {
    return(list(score = score, effect = "weak"))
  } else {
    return(list(score = score, effect = "neut"))
  }
}

reverseComplementMotif <- function(pwm) {
  rows <- rownames(pwm)
  cols <- colnames(pwm)
  Ns <- pwm["N", ]
  pwm <- pwm[4:1, length(cols):1]
  pwm <- rbind(pwm, Ns)
  rownames(pwm) <- rows
  colnames(pwm) <- cols
  return(pwm)
}


scoreSeqWindows <- function(ppm, seq) { #, collapse = TRUE) {
  ppm.width <- ncol(ppm)
  seq.len <- length(seq)
  diag.ind <- rep.int(ppm.width, seq.len - ppm.width)
  ranges <- vapply(c(0L, cumsum(diag.ind)),
                   function(x,
                            range = (1L + 0L:(ppm.width - 1L) * (ppm.width + 1L)))
                   {
                     x + range
                   },
                   integer(ppm.width))
  scores <- t(ppm[seq, ])[ranges]
  # scores_rc <- t(reverseComplementMotif(ppm)[seq, ])[ranges[, ncol(ranges):1]]
  scores_rc <- t(reverseComplementMotif(ppm)[seq, ])[ranges]
  scores <- split(scores, ceiling(seq_along(scores)/ppm.width))
  scores_rc <- split(scores_rc, ceiling(seq_along(scores_rc)/ppm.width))
  # if (collapse) {
  res <- vapply(Map(function(x, y) {matrix(data = c(x, y), nrow = 2,
                                           byrow = TRUE, dimnames = list(c(1, 2)))},
                    scores, scores_rc),
                rowSums,
                numeric(2))
  # } else {
  #   res <- Map(function(x, y) {matrix(data = c(x, y), nrow = 2,
  #                                     byrow = TRUE, dimnames = list(c(1, 2)))},
  #              scores, scores_rc)
  # }
  return(res)
}

maxThresholdWindows <- function(window.frame) {
  start.ind <- as.integer(colnames(window.frame)[1]) - 1L
  max.win <- arrayInd(which.max(window.frame), dim(window.frame))
  return(list(window = as.integer(colnames(window.frame)[max.win[, 2] + start.ind]),
              strand = c(1, 2)[max.win[, 1]]))
}


#' @import methods
#' @import GenomicRanges
#' @import S4Vectors
#' @import BiocGenerics
#' @import IRanges
#' @importFrom Biostrings getSeq replaceLetterAt reverseComplement complement replaceAt pairwiseAlignment insertion deletion matchPattern
#' @importFrom TFMPvalue TFMpv2sc
#' @importFrom stringr str_locate_all str_sub
scoreSnpList <- function(fsnplist, pwmList, method = "default", bkg = NULL,
                         threshold = 1e-3, show.neutral = FALSE, verbose = FALSE,
                         genome.bsgenome=NULL, pwmList.pc = NULL, pwmRanges = NULL, filterp=TRUE) {
  k <- max(sapply(pwmList, ncol))
  ref_len <- nchar(fsnplist$REF)
  alt_len <- nchar(fsnplist$ALT)
  is.indel <- ref_len > 1L | alt_len > 1L
  ## check that reference matches ref genome
  equals.ref <- getSeq(genome.bsgenome, fsnplist) == fsnplist$REF
  if (!all(equals.ref)) {
    stop(paste(names(fsnplist[!equals.ref]), "reference allele does not match value in reference genome",
               sep = " "))
  }
  if (sum(is.indel) < length(is.indel) & sum(is.indel) > 0L) {
    fsnplist.indel <- fsnplist[is.indel]
    fsnplist.snv <- fsnplist[!is.indel]
  } else if (sum(is.indel) == length(is.indel)) {
    fsnplist.indel <- fsnplist
    fsnplist.snv <- NULL
  } else if (sum(is.indel) == 0L) {
    fsnplist.indel <- NULL
    fsnplist.snv <- fsnplist
  }
  if (!is.null(fsnplist.indel)) {
    snp.sequence.ref.indel <- getSeq(genome.bsgenome, promoters(fsnplist.indel, upstream = k - 1,
                                                                downstream = k + max(nchar(fsnplist.indel$REF))))
    at <- as(IRanges(start = k, width = width(fsnplist.indel)), "IRangesList")
    snp.sequence.alt.indel <- DNAStringSet(Map(replaceAt,
                                               x = snp.sequence.ref.indel,
                                               at = at,
                                               fsnplist.indel$ALT))

    need.alignment <- !(lengths(fsnplist.indel$REF) == 1 | lengths(fsnplist.indel$ALT) == 1)
    insertion.var <- lengths(fsnplist.indel$REF) < lengths(fsnplist.indel$ALT)
    fsnplist.indel$ALT_loc <- 1L
    fsnplist.indel[insertion.var]$ALT_loc <- Map(seq,
                                           from = nchar(fsnplist.indel[insertion.var]$REF) + 1L,
                                           to = nchar(fsnplist.indel[insertion.var]$ALT))
    fsnplist.indel[!insertion.var]$ALT_loc <- Map(seq,
                                            from = nchar(fsnplist.indel[!insertion.var]$ALT) + 1L,
                                            to = nchar(fsnplist.indel[!insertion.var]$REF))
    ref.len <- nchar(fsnplist.indel[nchar(fsnplist.indel$REF) == nchar(fsnplist.indel$ALT)]$REF)
    fsnplist.indel[nchar(fsnplist.indel$REF) == nchar(fsnplist.indel$ALT)]$ALT_loc <- lapply(ref.len,
                                                                           function(x) {
                                                                             seq(from = 1L,
                                                                                 to = x)
                                                                           })
    fsnplist.indel$varType <- "Other"
    if (sum(insertion.var) > 0) {
      fsnplist.indel[insertion.var, ]$varType <- "Insertion"
    }
    if (sum(lengths(fsnplist.indel$REF) > lengths(fsnplist.indel$ALT)) > 0) {
      fsnplist.indel[lengths(fsnplist.indel$REF) > lengths(fsnplist.indel$ALT)]$varType <- "Deletion"
    }
    if (any(need.alignment)) {
      need.del <- any(!insertion.var & need.alignment)
      need.ins <- any(insertion.var & need.alignment)
      if (need.del) {
        pattern.del <- Map(matchPattern,
                           fsnplist.indel[!insertion.var & need.alignment]$ALT,
                           fsnplist.indel[!insertion.var & need.alignment]$REF,
                           with.indels = F, max.mismatch = 0)
        names(pattern.del) <- names(fsnplist.indel[!insertion.var & need.alignment])
        pattern.del <- sapply(pattern.del, function(x) {
          slen <- length(subject(x))
          x <- x[start(x) == 1 | end(x) == slen]
          if (length(x) > 0) {
            x <- x[1]
            x <- gaps(x)
            x <- as(x, "IRanges")
          }
        })
        pattern.del.valid <- sapply(pattern.del, function(x) {length(x) > 0})
        nr.pattern.del <- lapply(pattern.del[pattern.del.valid], function(x) {start(x):end(x)})
        fsnplist.indel[!insertion.var & need.alignment][names(pattern.del[pattern.del.valid])]$ALT_loc <- nr.pattern.del
        if (any((!insertion.var & need.alignment)[!pattern.del.valid])) {
          alignment.del <- Map(pairwiseAlignment,
                               fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid]$ALT,
                               fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid]$REF,
                               type = "global")
          names(alignment.del) <- names(fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid])
          alignment.del <- sapply(sapply(alignment.del, deletion), unlist)
          alignment.del.valid <- sapply(alignment.del, function(x) {length(x) > 0})
          nr.alignment.del <- lapply(alignment.del[alignment.del.valid], function(x) {start(x):end(x)})
          fsnplist.indel[!insertion.var & need.alignment][names(alignment.del[alignment.del.valid])]$ALT_loc <- nr.alignment.del
          rm(alignment.del, nr.alignment.del)
        }
        rm(pattern.del, nr.pattern.del, pattern.del.valid, need.del)
      }
      if (need.ins) {
        pattern.ins <- Map(matchPattern,
                           fsnplist.indel[insertion.var & need.alignment]$REF,
                           fsnplist.indel[insertion.var & need.alignment]$ALT,
                           with.indels = F, max.mismatch = 0)
        names(pattern.ins) <- names(fsnplist.indel[insertion.var & need.alignment])
        pattern.ins <- sapply(pattern.ins, function(x) {
          slen <- length(subject(x))
          x <- x[start(x) == 1 | end(x) == slen]
          if (length(x) > 0) {
            x <- x[1]
            x <- gaps(x)
            x <- as(x, "IRanges")
          }
        })
        pattern.ins.valid <- sapply(pattern.ins, function(x) {length(x) > 0})
        nr.pattern.ins <- lapply(pattern.ins[pattern.ins.valid], function(x) {start(x):end(x)})
        fsnplist.indel[insertion.var & need.alignment][names(pattern.ins[pattern.ins.valid])]$ALT_loc <- nr.pattern.ins
        if (any((insertion.var & need.alignment)[!pattern.ins.valid])) {
          alignment.ins <- Map(pairwiseAlignment,
                               fsnplist.indel[insertion.var & need.alignment][!pattern.ins.valid]$ALT,
                               fsnplist.indel[insertion.var & need.alignment][!pattern.ins.valid]$REF,
                               type = "global")
          names(alignment.ins) <- names(fsnplist.indel[insertion.var & need.alignment][!pattern.ins.valid])
          alignment.ins <- sapply(sapply(alignment.ins, insertion), unlist)
          alignment.ins.valid <- sapply(alignment.ins, function(x) {length(x) > 0})
          nr.alignment.ins <- lapply(alignment.ins[alignment.ins.valid], function(x) {start(x):end(x)})
          fsnplist.indel[insertion.var & need.alignment][names(alignment.ins[alignment.ins.valid])]$ALT_loc <- nr.alignment.ins
          rm(alignment.ins, nr.alignment.ins)
        }
        rm(pattern.ins, nr.pattern.ins, pattern.ins.valid, need.ins)
      }
      rm(ref.len, insertion.var, need.alignment)
    }
  }
  if (!is.null(fsnplist.snv)) {
    snp.sequence.ref.snv <- getSeq(genome.bsgenome, promoters(fsnplist.snv, upstream = k - 1,
                                                              downstream = k + 1))
    at <- matrix(FALSE, nrow = length(snp.sequence.ref.snv), ncol = (k * 2))
    at[, k] <- TRUE
    snp.sequence.alt.snv <- replaceLetterAt(snp.sequence.ref.snv, at, fsnplist.snv$ALT)
    fsnplist.snv$ALT_loc <- 1L
    fsnplist.snv$varType <- "SNV"
  }
  if (sum(is.indel) < length(is.indel) & sum(is.indel) > 0L) {
    fsnplist <- c(fsnplist.indel, fsnplist.snv)
    snp.sequence.alt <- strsplit(as.character(c(snp.sequence.alt.indel,
                                                snp.sequence.alt.snv)), "")
    snp.sequence.ref <- strsplit(as.character(c(snp.sequence.ref.indel,
                                                snp.sequence.ref.snv)), "")
    rm(fsnplist.indel, fsnplist.snv,
       snp.sequence.alt.indel, snp.sequence.ref.indel,
       snp.sequence.alt.snv, snp.sequence.ref.snv)
  } else if (sum(is.indel) == length(is.indel)) {
    fsnplist <- fsnplist.indel
    snp.sequence.alt <- strsplit(as.character(snp.sequence.alt.indel), "")
    snp.sequence.ref <- strsplit(as.character(snp.sequence.ref.indel), "")
    rm(fsnplist.indel, snp.sequence.alt.indel, snp.sequence.ref.indel)
  } else if (sum(is.indel) == 0L) {
    fsnplist <- fsnplist.snv
    snp.sequence.alt <- strsplit(as.character(snp.sequence.alt.snv), "")
    snp.sequence.ref <- strsplit(as.character(snp.sequence.ref.snv), "")
    rm(fsnplist.snv, snp.sequence.alt.snv, snp.sequence.ref.snv)
  }
  rm(at); gc()

  res.el.e <- new.env()
  for (snp.map.i in seq_along(snp.sequence.alt)) {
    snp.ref <- snp.sequence.ref[[snp.map.i]]
    snp.alt <- snp.sequence.alt[[snp.map.i]]
    ref.len <- nchar(fsnplist[snp.map.i]$REF)
    alt.len <- nchar(fsnplist[snp.map.i]$ALT)
    alt.loc <- fsnplist[snp.map.i]$ALT_loc[[1]]
    res.el <- rep(fsnplist[snp.map.i], length(pwmList))
    res.el$motifPos <- as.integer(NA)
    res.el$motifID <- mcols(pwmList)$providerID
    res.el$geneSymbol <- mcols(pwmList)$geneSymbol
    res.el$dataSource <- mcols(pwmList)$dataSource
    res.el$providerName <- mcols(pwmList)$providerName
    res.el$providerId <- mcols(pwmList)$providerId
    res.el$seqMatch <- as.character(NA)
    res.el$pctRef <- as.numeric(NA)
    res.el$pctAlt <- as.numeric(NA)
    if (filterp) {
      res.el$scoreRef <- as.numeric(NA)
      res.el$scoreAlt <- as.numeric(NA)
      res.el$Refpvalue <- as.numeric(NA)
      res.el$Altpvalue <- as.numeric(NA)
    }
    if (ref.len > 1 | alt.len > 1) {
      res.el$altPos <- as.numeric(NA)
      res.el$alleleDiff <- as.numeric(NA)
    } else {
      res.el$snpPos <- as.integer(NA)
      res.el$alleleRef <- as.numeric(NA)
      res.el$alleleAlt <- as.numeric(NA)
    }
    res.el$effect <- as.character(NA)
    for (pwm.i in seq_along(pwmList)) {
      pwm.basic <- pwmList[[pwm.i]]
      pwm <- pwmList.pc[[pwm.i]]
      len <- ncol(pwm)
      thresh <- threshold[[pwm.i]]
      if (ref.len > 1 | alt.len > 1) {
        seq.len <- max(alt.loc)
        seq.start <- min(alt.loc)
        alt.range <- ref.range <- (k - (ncol(pwm) - seq.start)):(k + ncol(pwm) + seq.len - 2)
        if (!show.neutral & identical(snp.ref[ref.range], snp.alt[alt.range])) next()
      } else {
        alt.range <- ref.range <- (k - (ncol(pwm) - 1)):(k + ncol(pwm) - 1)
      }
      seq.remove <- ref.len - alt.len
      if (seq.remove < 0) {
        ref.range <- ref.range[1:(length(ref.range) + seq.remove)]
      } else {
        alt.range <- alt.range[1:(length(alt.range) - seq.remove)]
      }
      ref.windows <- scoreSeqWindows(ppm = pwm, seq = snp.ref[ref.range])
      alt.windows <- scoreSeqWindows(ppm = pwm, seq = snp.alt[alt.range])
      if (!filterp) {
        ref.windows <- (ref.windows - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1])
        alt.windows <- (alt.windows - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1])
      }
      if (any(alt.windows > thresh) | any(ref.windows > thresh)) {
        hit.alt <- maxThresholdWindows(alt.windows)
        hit.ref <- maxThresholdWindows(ref.windows)
        bigger <- ref.windows[hit.ref$strand, hit.ref$window] >= alt.windows[hit.alt$strand, hit.alt$window]
        if (bigger) {
          hit <- hit.ref
        } else {
          hit <- hit.alt
        }
      } else {
        hit.alt <- list(window = 0L, strand = 0L)
        hit.ref <- list(window = 0L, strand = 0L)
        hit <- NULL
      }
      if (!show.neutral) {
        if (identical(alt.windows[hit.alt$strand, hit.alt$window],
                      ref.windows[hit.ref$strand, hit.ref$window])) next()
      }
      if (!is.null(hit)) {
        result <- res.el[pwm.i]
        uniquename <- paste(names(result), result$dataSource, result$providerName, result$providerId, sep = "%%")
        if (nchar(result$REF) > 1 | nchar(result$ALT) > 1) {
          # indel
          scorediff <- scoreIndel(pwm = pwm,
                                  ref.seq = snp.ref[ref.range],
                                  alt.seq = snp.alt[alt.range],
                                  hit.ref = hit.ref, hit.alt = hit.alt)
          effect <- scorediff$effect
          score <- scorediff$score
          #ref.pos <- ref.range[(ncol(pwm) - (seq.start - 1)):((ncol(pwm) + ref.len - seq.start))]
          ref.pos <- k:(k + nchar(result$REF) - 1L)
          #alt.pos <- alt.range[(ncol(pwm) - (seq.start - 1)):((ncol(pwm) + alt.len - seq.start))]
          alt.pos <- k:(k + nchar(result$ALT) - 1L)
          if (effect == "neut") {
            if (show.neutral) {
              res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                           snp.ref, snp.alt,
                                                           ref.pos, alt.pos,
                                                           hit.ref, hit.alt,
                                                           ref.windows, alt.windows,
                                                           score, effect, len,
                                                           k, pwm, calcp = filterp)
            }
          } else {
            res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                         snp.ref, snp.alt,
                                                         ref.pos, alt.pos,
                                                         hit.ref, hit.alt,
                                                         ref.windows, alt.windows,
                                                         score, effect, len,
                                                         k, pwm, calcp = filterp)
          }
        } else {
          snp.pos <- len - hit$window + 1L
          if (hit$strand == 1) {
            allelR <- pwm.basic[as.character(result$REF), snp.pos]
            allelA <- pwm.basic[as.character(result$ALT), snp.pos]
          } else {
            allelR <- pwm.basic[as.character(complement(result$REF)), snp.pos]
            allelA <- pwm.basic[as.character(complement(result$ALT)), snp.pos]
          }
          scorediff <- snpEff(allelR, allelA)
          effect <- scorediff$effect
          score <- scorediff$score
          if (effect == "neut") {
            if (show.neutral) {
              res.el.e[[uniquename]] <- updateResultsSnv(result, snp.ref[ref.range], snp.pos,
                                                         hit, ref.windows, alt.windows,
                                                         allelR, allelA, effect, len,
                                                         k, pwm, calcp = filterp)
            }
          } else {
            res.el.e[[uniquename]] <- updateResultsSnv(result, snp.ref[ref.range], snp.pos,
                                                       hit, ref.windows, alt.windows,
                                                       allelR, allelA, effect, len,
                                                       k, pwm, calcp = filterp)
          }
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
    if ("ALT_loc" %in% names(mcols(resultSet))) mcols(resultSet)$ALT_loc <- NULL
    max.match <- max(vapply(str_locate_all(resultSet$seqMatch, "\\w"), max, integer(1)))
    min.match <- min(vapply(str_locate_all(resultSet$seqMatch, "\\w"), min, integer(1)))
    resultSet$seqMatch <- str_sub(resultSet$seqMatch,
                                  start = min.match + 1,
                                  end = max.match + 1)
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
updateResultsSnv <- function(result, snp.seq, snp.pos, hit, ref.windows, alt.windows,
                          allelR, allelA, effect, len, k, pwm, calcp) {
  strand.opt <- c("+", "-")
  strand(result) <- strand.opt[[hit$strand]]
  hit$window <- as.integer(hit$window)
  mresult <- mcols(result)
  mresult[["snpPos"]] <- start(result)
  mresult[["motifPos"]] <- as.integer(snp.pos)
  matchs <- snp.seq
  seq.pos <- snp.pos + hit$window - 1
  # matchs <- snp.seq[(k - len + (hit$window)):(k + (hit$window) - 1)]
  matchs[-(seq.pos)] <- tolower(matchs[-(seq.pos)])
  matchs <- paste(matchs, collapse = "")
  mresult[["seqMatch"]] <- str_pad(matchs, width = k * 2, side = "both")
  start(result) <- start(result) - snp.pos + 1
  end(result) <- end(result) - snp.pos + len
  if (calcp) {
    mresult[["scoreRef"]] <- ref.windows[hit$strand, hit$window]
    mresult[["scoreAlt"]] <- alt.windows[hit$strand, hit$window]
    mresult[["Refpvalue"]] <- NA
    mresult[["Altpvalue"]] <- NA
    pwmrange <- colSums(colRanges(pwm))
    mresult[["pctRef"]] <- (mresult[["scoreRef"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
    mresult[["pctAlt"]] <- (mresult[["scoreAlt"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
  } else {
    mresult[["pctRef"]] <- ref.windows[hit$strand, hit$window]
    mresult[["pctAlt"]] <- alt.windows[hit$strand, hit$window]
  }
  mresult[["alleleRef"]] <- allelR
  mresult[["alleleAlt"]] <- allelA
  mresult[["effect"]] <- effect
  mcols(result) <- mresult
  return(result)
}


updateResultsIndel <- function(result,
                               ref.seq, alt.seq,
                               ref.pos, alt.pos,
                               hit.ref, hit.alt,
                               ref.windows, alt.windows,
                               score, effect, len, k, pwm, calcp) {
  strand.opt <- c("+", "-")
  if (score > 0L) {
    best.hit <- hit.alt
    matchs <- alt.seq
    snp.pos <- alt.pos
  } else {
    best.hit <- hit.ref
    matchs <- ref.seq
    snp.pos <- ref.pos
  }
  strand(result) <- strand.opt[[best.hit$strand]]
  best.hit$window <- as.integer(best.hit$window)
  mresult <- mcols(result)
  # if (mresult$SNP_id == "rs572297057;rs541626760" & mresult$providerName %in% c("MA0070.1", "MA1114.1", "MA1113.1")) browser()
  alt_loc <- range(mresult$ALT_loc)
  ref_start <- (1 - alt_loc[[1]])
  ref_start <- ifelse(ref_start <= 0, ref_start - 1, ref_start)
  # motif.start <- (alt_loc[[1]] - 1) + (-len) + (best.hit$window - 1)
  motif.start <- (alt_loc[[1]]) + (-len) + (best.hit$window) + ref_start
  motif.start <- ifelse(motif.start >= 0, motif.start + 1, motif.start)
  if ((mresult$varType == "Insertion" & score < 0) |
      (mresult$varType == "Deletion" & score > 0)) {
    motif.end <- motif.start + len
  } else {
    if (motif.start > 0) {
      motif.end <- len - length(motif.start:length(alt_loc[1]:alt_loc[2]))
    } else {
      motif.end <- motif.start + len - length(alt_loc[1]:alt_loc[2])
    }
  }
  motif.end <- ifelse(motif.end <= 0, motif.end - 1, motif.end)
  mresult$motifPos <- list(c(motif.start, motif.end))
  mresult$altPos <- mresult$ALT_loc
  seq.range <- (k - (len - alt_loc[[1]])):(k + len + alt_loc[[2]] - 2)
  matchs[-(snp.pos)] <- tolower(matchs[-(snp.pos)])
  matchs <- paste(matchs[seq.range], collapse = "")
  mresult[["seqMatch"]] <- str_pad(matchs, width = (k * 2) + alt_loc[[2]], side = "both")
  if (calcp) {
    mresult[["scoreRef"]] <- ref.windows[hit.ref$strand, hit.ref$window]
    mresult[["scoreAlt"]] <- alt.windows[hit.alt$strand, hit.alt$window]
    mresult[["Refpvalue"]] <- NA
    mresult[["Altpvalue"]] <- NA
    pwmrange <- colSums(colRanges(pwm[-5,]))
    mresult[["pctRef"]] <- (mresult[["scoreRef"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
    mresult[["pctAlt"]] <- (mresult[["scoreAlt"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
  } else {
    mresult[["pctRef"]] <- ref.windows[hit.ref$strand, hit.ref$window]
    mresult[["pctAlt"]] <- alt.windows[hit.alt$strand, hit.alt$window]
  }
  mresult[["alleleDiff"]] <- score
  mresult[["effect"]] <- effect
  mcols(result) <- mresult
  return(result)
}

preparePWM <- function(pwmList = pwmList,
                       filterp = filterp,
                       bkg = bkg,
                       scoreThresh = threshold,
                       method = "default") {

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
    pwmThresh <- lapply(pwmList.pc2, TFMpv2sc, pvalue = scoreThresh, bg = bkg, type = "PWM")
    pwmThresh <- Map("+", pwmThresh, -0.02)
  } else {
    pwmRanges <- Map(function(pwm, omega) {
      x <- colSums(colRanges(pwm))
      return(x)
    }, pwmList.pc, pwmOmegas)
    pwmThresh <- rep.int(scoreThresh, times = length(pwmRanges))
  }
  pwmList@listData <- lapply(pwmList, function(pwm) { pwm <- rbind(pwm, N=0); return(pwm) })
  pwmList.pc <- lapply(pwmList.pc, function(pwm) { pwm <- rbind(pwm, N=0); return(pwm) })
  return(list(pwmList = pwmList,
              pwmListPseudoCount = pwmList.pc,
              pwmRange = pwmRanges,
              pwmThreshold = pwmThresh))
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
#' @importFrom BiocParallel bplapply
#' @importFrom stringr str_length str_trim
#' @export
motifbreakR <- function(snpList, pwmList, threshold=0.85, filterp = FALSE,
                        method = "default", show.neutral = FALSE, verbose = FALSE,
                        bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                        BPPARAM=bpparam()) {
  ## Cluster / MC setup
  if (.Platform$OS.type == "windows" && inherits(BPPARAM, "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n",
            "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpnworkers(BPPARAM)
  num.snps <- length(snpList)
  if (num.snps < cores) {
    cores <- num.snps
  }
  if (is(BPPARAM, "SnowParam")) {
    bpstart(BPPARAM)
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  ##
  ## Genome Setup
  genome.package <- attributes(snpList)$genome.package
  if (requireNamespace(eval(genome.package), quietly = TRUE, character.only = TRUE)) {
    genome.bsgenome <- eval(parse(text = paste(genome.package, genome.package, sep="::")))
  } else {
    stop(paste0(eval(genome.package), " is the genome selected for this snp list and \n",
                "is not present on your environment. Please load it and try again."))
  }
  ##

  snpList <- sapply(suppressWarnings(split(snpList, 1:cores)), list)

  pwms <- preparePWM(pwmList = pwmList, filterp = filterp,
                     scoreThresh = threshold, bkg = bkg,
                     method = method)

  # x <- lapply(snpList, scoreSnpList,
  #             pwmList = pwms$pwmList, threshold = pwms$pwmThreshold,
  #             pwmList.pc = pwms$pwmListPseudoCount, pwmRanges = pwms$pwmRange,
  #             method = method, bkg = bkg, show.neutral = show.neutral,
  #             verbose = ifelse(cores == 1, verbose, FALSE), genome.bsgenome = genome.bsgenome,
  #             filterp = filterp)
  x <- try(bplapply(snpList, scoreSnpList,
                    pwmList = pwms$pwmList, threshold = pwms$pwmThreshold,
                    pwmList.pc = pwms$pwmListPseudoCount, pwmRanges = pwms$pwmRange,
                    method = method, bkg = bkg, show.neutral = show.neutral,
                    verbose = ifelse(cores == 1, verbose, FALSE), genome.bsgenome = genome.bsgenome,
                    filterp = filterp, BPPARAM = BPPARAM))

  ## Cluster / MC cleanup
  if (inherits(x, "try-error")) {
    if (is(BPPARAM, "SnowParam")) {
      bpstop(BPPARAM)
    }
    stop(attributes(x)$condition)
  }
  if (is(BPPARAM, "SnowParam")) {
	  bpstop(BPPARAM)
  }

  drops <- sapply(x, is.null)
  x <- x[!drops]
  pwmList <- pwms$pwmList
  pwmList@listData <- lapply(pwms$pwmList, function(pwm) { pwm <- pwm[c("A", "C", "G", "T"), ]; return(pwm) })
  pwmList.pc <- lapply(pwms$pwmListPseudoCount, function(pwm) { pwm <- pwm[c("A", "C", "G", "T"), ]; return(pwm) })

  if (length(x) > 1) {
    x <- unlist(GRangesList(unname(x)))
    snpList <- unlist(GRangesList(snpList), use.names = F)
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
    pvalues <- lapply(seq_along(results), function(i, pwmList, pwmListmeta, bkg) {
                        result <- results[i]
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


# #' @importFrom grImport PostScriptTrace readPicture grid.picture
# addPWM.stack <- function(identifier, index, GdObject, ...) {
#   psloc <- tempdir()
#   ps.file <- paste(psloc, "stack.ps", sep = "/")
#   PostScriptTrace(ps.file, paste0(ps.file, ".xml"))
#   motif.figure <- readPicture(paste0(ps.file, ".xml"))
#   # motif.figure <- motif.figure[-1]
#   snppos <- sapply(sapply(mcols(GdObject@range)[, "feature"], strsplit, "@"), "[[", 2)
#   motif.i <- 1
#   highlight <- snppos[[motif.i]]
#   paths <- motif.figure@paths
#   for(path.i in seq_along(paths)) {
#     if (inherits(paths[[path.i]], "PictureText")) {
#       if (paths[[path.i]]@string == highlight) {
#         path <- paths[[path.i]]
#         path@rgb <- "#FF0000"
#         path@letters <- lapply(path@letters, function(letters) {
#           letters@rgb <- "#FF0000"
#           return(letters)
#         })
#         paths[[path.i]] <- path
#         motif.i <- motif.i + 1
#         if(motif.i <= length(snppos)) {
#           highlight <- snppos[[motif.i]]
#         } else {
#           motif.i <- length(snppos)
#         }
#       }
#     }
#   }
#   motif.figure@paths <- paths
#   grid.picture(motif.figure, distort = FALSE)
# }

#' @importFrom grImport2 readPicture grid.picture
addPWM.stack <- function(identifier, index, GdObject, pwm_set, ...) {
  # motifStack(pwm_set, ncex = 1.0, layout = "stack")
  psloc <- tempdir()
  # ps.file <- paste(psloc, "stack.ps", sep = "/")
  # PostScriptTrace(ps.file, paste0(ps.file, ".xml"))
  ps.file <- paste(psloc, "stack.svg", sep = "/")
  motif.figure <- grImport2::readPicture(ps.file)
  # # motif.figure <- motif.figure[-1]
  # snppos <- sapply(sapply(mcols(GdObject@range)[, "feature"], strsplit, "@"), "[[", 2)
  # motif.i <- 1
  # highlight <- snppos[[motif.i]]
  # paths <- motif.figure@paths
  # for(path.i in seq_along(paths)) {
  #   if (inherits(paths[[path.i]], "PictureText")) {
  #     if (paths[[path.i]]@string == highlight) {
  #       path <- paths[[path.i]]
  #       path@rgb <- "#FF0000"
  #       path@letters <- lapply(path@letters, function(letters) {
  #         letters@rgb <- "#FF0000"
  #         return(letters)
  #       })
  #       paths[[path.i]] <- path
  #       motif.i <- motif.i + 1
  #       if(motif.i <= length(snppos)) {
  #         highlight <- snppos[[motif.i]]
  #       } else {
  #         motif.i <- length(snppos)
  #       }
  #     }
  #   }
  # }
  # motif.figure@paths <- paths
  grImport2::grid.picture(motif.figure, distort = FALSE)
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

plotMotifLogoStack.2 <- function(pfms, ...) {
  pfms <- rev(pfms)
  n <- length(pfms)
  lapply(pfms, function(.ele) {
    if (!is(.ele, "pfm"))
      stop("pfms must be a list of class pfm")
  })
  opar <- par(mfrow = c(n, 1), mar = c(3.5, 3.5, 1.5, 0.5))
  assign("tmp_motifStack_symbolsCache", list(), pos = ".GlobalEnv")
  motifStack::plotMotifLogo(pfms[[1]], motifName = pfms[[1]]@name, p = rep(0.25, 4))
  for (i in seq.int(n)[-1]) {
    motifStack::plotMotifLogo(pfms[[n]], motifName = pfms[[n]]@name,
                              p=rep(0.25, 4), xlab = NA, newpage = FALSE)
  }
  rm(list = "tmp_motifStack_symbolsCache", pos = ".GlobalEnv")
  par(opar)
}

#' @importFrom grid grid.newpage pushViewport viewport popViewport
plotMotifLogoStack.3 <- function(pfms, ...) {
  n <- length(pfms)
  lapply(pfms, function(.ele) {
    if (class(.ele) != "pfm")
      stop("pfms must be a list of class pfm")
  })
  assign("tmp_motifStack_symbolsCache", list(), pos = ".GlobalEnv")
  grid.newpage()
  ht <- 1/n
  y0 <- 0.5 * ht
  for (i in rev(seq.int(n))) {
    pushViewport(viewport(y = y0, height = ht))
    plotMotifLogo(pfms[[i]], motifName = pfms[[i]]@name, ncex = 1,
                  p = pfms[[i]]@background, colset = pfms[[i]]@color,
                  xlab = NA, newpage = FALSE, margins = c(1.5, 4.1,
                                                          1.1, 0.1), ...)
    popViewport()
    y0 <- y0 + ht
  }
  rm(list = "tmp_motifStack_symbolsCache", pos = ".GlobalEnv")
  return()
}

#' @importFrom stringr str_replace
#' @importFrom motifStack addBlank
DNAmotifAlignment.2snp <- function(pwms, result) {
  from <- min(sapply(result$motifPos, `[`, 1))
  to <- max(sapply(result$motifPos, `[`, 2))
#  pos <- mcols(result)$motifPos
#  pos <- pos[as.logical(strand(result) == "+")][1]
  for (pwm.i in seq_along(pwms)) {
    #pwm <- pwms[[pwm.i]]@mat
    ## get pwm info from result data
    pwm.name <- pwms[[pwm.i]]@name
    pwm.name <- str_replace(pwm.name, pattern = "-:rc$", replacement = "")
    pwm.name <- str_replace(pwm.name, pattern = "-:r$", replacement = "")
    pwm.info <- attributes(result)$motifs
    pwm.id <- mcols(pwm.info[pwm.name, ])$providerId
    pwm.name <- mcols(pwm.info[pwm.name, ])$providerName
    mresult <- result[result$providerId == pwm.id & result$providerName == pwm.name, ]
    mstart <- mresult$motifPos[[1]][1]
    mend <- mresult$motifPos[[1]][2]
    # browser()
    if ((mcols(mresult)$varType == "Insertion" & mcols(mresult)$alleleDiff < 0) |
        (mcols(mresult)$varType == "Deletion" & mcols(mresult)$alleleDiff > 0)) {
      new.mat <- cbind(pwms[[pwm.i]]@mat[, 1:abs(mresult$motifPos[[1]][1])],
                       matrix(c(0.25, 0.25, 0.25, 0.25), ncol = length(mcols(mresult)$altPos[[1]]), nrow = 4),
                       pwms[[pwm.i]]@mat[, (abs(mresult$motifPos[[1]][1]) + 1):ncol(pwms[[pwm.i]]@mat)])
      pwms[[pwm.i]]@mat <- new.mat
      start.offset <- mstart - from
      end.offset <- to - mend
    } else {
      if (mstart < 0 | from > 0) {
        start.offset <- mstart - from
      } else {
        start.offset <- (mstart - 1) - from
      }
      if (mend > 0 | to < 0) {
        end.offset <- to - mend
      } else {
        end.offset <- to - (mend + 1)
      }
    }
    if (start.offset > 0) {
      pwms[[pwm.i]] <- addBlank(x = pwms[[pwm.i]], n = start.offset, b = FALSE)
    }
    if (end.offset > 0) {
      pwms[[pwm.i]] <- addBlank(x = pwms[[pwm.i]], n = end.offset, b = TRUE)
    }
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
  result <- result[order(sapply(result$motifPos, min), sapply(result$motifPos, max)), ]
  result <- result[result$effect %in% effect]
  chromosome <- as.character(seqnames(result))[[1]]
  genome.package <- attributes(result)$genome.package
  genome.bsgenome <- eval(parse(text = genome.package))
  seq.len <- max(length(result$REF[[1]]), length(result$ALT[[1]]))
  distance.to.edge <- max(abs(c(sapply(result$motifPos, min),
                                sapply(result$motifPos, max)))) + 4
  from <- start(result)[[1]] - distance.to.edge + 1
  to <- end(result)[[1]] + distance.to.edge
  pwmList <- attributes(result)$motifs
  pwm.names <- result$providerId
  if (stackmotif) {
    getmotifs <- mcols(pwmList)$providerId %in% result$providerId &
      mcols(pwmList)$providerName %in% result$providerName
    pwms <- pwmList[getmotifs, ]
    pwms <- pwms[order(match(paste0(mcols(pwms)$providerId, mcols(pwms)$providerName),
                             paste0(result$providerId, result$providerName)))]
    if (reverseMotif) {
      for (pwm.i in seq_along(pwms)) {
        pwm.name <- names(pwms[pwm.i])
        pwm.id <- mcols(pwms[pwm.name, ])$providerId
        pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
        doRev <- as.logical(strand(result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]) == "-")
        if (doRev) {
          pwm <- pwms[[pwm.i]]
          pwm <- pwm[, rev(1:ncol(pwm))]
          rownames(pwm) <- c("T", "G", "C", "A")
          pwm <- pwm[c("A", "C", "G", "T"), ]
          pwms[[pwm.i]] <- pwm
          names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], "-:rc")
        }
      }
    }
    # else {
    #   for (pwm.i in seq_along(pwms)) {
    #     pwm.name <- names(pwms[pwm.i])
    #     pwm.id <- mcols(pwms[pwm.name, ])$providerId
    #     pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
    #     doRev <- as.logical(strand(result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]) == "-")
    #     if(doRev) {
    #       pwm <- pwms[[pwm.i]]
    #       pwm <- pwm[, rev(1:ncol(pwm))]
    #       pwms[[pwm.i]] <- pwm
    #       names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], "-:r")
    #     }
    #   }
    # }
    pwms <- lapply(names(pwms), function(x, pwms=pwms) {new("pfm", mat = pwms[[x]],
                                                            name = x)}, pwms)
    pwms <- DNAmotifAlignment.2snp(pwms, result)
    pwmwide <- max(sapply(pwms, function(x) { ncol(x@mat)}))

    psloc <- tempdir()
    old.list <- dev.list()
    svg(paste(psloc, "stack.svg", sep = "/"),
        width = pwmwide * (1/3),
        height = 2 * length(pwm.names))
    new.devlist <- dev.list()[!dev.list() %in% old.list]
    Sys.sleep(1)
    dev.set(new.devlist)
    # postscript(paste(psloc, "stack.ps", sep = "/"),
    #            width = pwmwide * (1/3),
    #            height = 2 * length(pwm.names),
    #            paper = "special", horizontal = FALSE,
    #            fonts = c("sans"))
    # markerStarts <- sapply(result$motifPos, `[`, 1)
    # markerEnds <- sapply(result$motifPos, `[`, 2)
    # markerRect <- new("marker", type = "rect",
    #                   start = markerStarts,
    #                   stop = markerEnds,
    #                   gp = gpar(lty = 2, fill = NA, col = "red"))
    # sapply(pwms, function(x) {ncol(x@mat)})
    theplot <- try(plotMotifLogoStack.3(pwms))
    if (inherits(theplot, "try-error")) {
      dev.off()
      stop("error plotting motifs: use graphics.off() prior to plotting")
    }
    dev.off()
    if (!is.null(dev.list())) dev.set(max(dev.list()))
  }
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
  if (result$varType[[1]] == "Deletion") {
    reflen <- length(result$REF[[1]])
    addedN <- DNAString(paste0(rep.int(".", reflen), collapse = ""))
    addedN <- replaceLetterAt(addedN, at = (1:reflen)[-result$altPos[[1]]], result$ALT[[1]])
    axisT <- GenomeAxisTrack(exponent = 0)
    seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", "auto"))
    altseq <- replaceAt(x = altseq, at = at, addedN)
  } else if (result$varType[[1]] == "Insertion") {
    # browser()
    altlen <- length(result$ALT[[1]])
    addedN <- DNAString(paste0(rep.int(".", altlen), collapse = ""))
    addedN <- replaceLetterAt(addedN, at = (1:altlen)[-result$altPos[[1]]], result$REF[[1]])
    refseq <- genome.bsgenome[[chromosome]]
    refseq <- DNAStringSet(replaceAt(x = refseq, at = at, addedN))
    altseq <- replaceAt(x = altseq, at = at, result$ALT[[1]])
    names(refseq) <- chromosome
    seqT <- SequenceTrack(refseq,
                          fontcolor = c(colorset("DNA", "auto"), N = "#FFFFFF", . = "#FFE3E6"),
                          chromosome = chromosome)
  } else {
    axisT <- GenomeAxisTrack(exponent = 0)
    altseq <- replaceAt(x = altseq, at = at, result$ALT[[1]])
    seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", "auto"))
  }
  altseq <- DNAStringSet(altseq)
  names(altseq) <- chromosome
  seqAltT <- SequenceTrack(altseq,
                           fontcolor = c(colorset("DNA", "auto"), N = "#FFFFFF", . = "#FFE3E6"),
                           chromosome = chromosome)


  ### end replace longer sections


  #altseq <- replaceLetterAt(altseq, at = wherereplace, letter = rep.int("N", sum(wherereplace)))
  #altseq <- replaceLetterAt(altseq, at = !wherereplace, letter = result$ALT[[1]])
  hiT <- HighlightTrack(trackList = list(seqT, seqAltT),
                        start = start(result[1]) + min(result[1]$altPos[[1]]) - 1,
                        end = start(result[1]) + min(result[1]$altPos[[1]]) - 2 +
                          length(result[1]$altPos[[1]]),
                        chromosome = chromosome)
  if (stackmotif) {
    selectingfun <- selcor
    detailfun <- addPWM.stack
  } else {
    detailfun <- addPWM.stack
    selectingfun <- selall
  }

  getmotifs <- mcols(pwmList)$providerId %in% result$providerId & mcols(pwmList)$providerName %in% result$providerName
  motif_ids <- names(pwmList)[getmotifs]
  names(motif_ids) <- mcols(pwmList)$providerName[getmotifs]

  for (mymotif_i in seq_along(result)) {
    mymotif <- result[mymotif_i]
    start(mymotif) <- start(mymotif) + min(mymotif$altPos[[1]]) - 1
    width(mymotif) <- length(mymotif$altPos[[1]])
    variant.start <- start(mymotif)
    variant.end <- end(mymotif)
    if (mymotif$motifPos[[1]][1] < 0) {
      start(mymotif) <- start(mymotif) + (mymotif$motifPos[[1]][1])
    } else {
      start(mymotif) <- start(mymotif) + (mymotif$motifPos[[1]][1] - 1)
    }
    if (mymotif$motifPos[[1]][2] < 0) {
      end(mymotif) <- end(mymotif) + (mymotif$motifPos[[1]][2] + 1)
    } else {
      end(mymotif) <- end(mymotif) + (mymotif$motifPos[[1]][2])
    }
    if ((result[mymotif_i]$varType == "Deletion" & result[mymotif_i]$alleleDiff > 0) |
        (result[mymotif_i]$varType == "Insertion" & result[mymotif_i]$alleleDiff < 0)) {
      mymotif <- c(mymotif, mymotif)
      end(mymotif)[1] <- variant.start - 1
      start(mymotif)[2] <- variant.end + 1
      mymotif[which.min(width(mymotif))]$motifPos <- NA
    }
    if (exists("mres")) {
      mres <- c(mres, mymotif)
    } else {
      mres <- mymotif
    }
  }
  result <- mres; rm(mres)
  motif_ids <- motif_ids[result$providerName]
  motifT <- AnnotationTrack(result, id = motif_ids,
                            fun = detailfun, group = result$providerName,
                            # feature = paste0("snp@",
                            #                 result$motifPos) ,
                            feature = ifelse(!is.na(result$motifPos), paste0(result$geneSymbol, "_motif"), ""),
                            name = names(result)[[1]], selectFun = selectingfun)
  if (exists("axisT")) {
    track_list <- list(ideoT, motifT, hiT, axisT)
  } else {
    track_list <- list(ideoT, motifT, hiT)
  }
  plotTracks(track_list, from = from, to = to, showBandId = TRUE,
             cex.main = 0.8, col.main = "darkgrey",
             add53 = TRUE, labelpos = "below", chromosome = chromosome, groupAnnotation = "group",
             collapse = FALSE, min.width = 1, featureAnnotation = "feature", cex.feature = 0.8,
             details.size = ifelse(stackmotif, 0.85, 0.5), detailsConnector.pch = NA,
             detailsConnector.lty = ifelse(stackmotif, 0, 3),
             shape = "box", cex.group = 0.8, fonts = c("sans", "Helvetica"))
  return(invisible(NULL))
}

