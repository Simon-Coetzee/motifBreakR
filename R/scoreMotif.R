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

prepareVariants <- function(fsnplist, genome.bsgenome, max.pwm.width) {
  k <- max.pwm.width; rm(max.pwm.width)
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
                           with.indels = FALSE, max.mismatch = 0)
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
          alt.width <- width(fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid]$ALT)
          ref.width <- width(fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid]$REF)
          names(alignment.del) <- names(fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid])
          alignment.ins <- sapply(sapply(alignment.del, insertion), unlist)
          alignment.del <- sapply(sapply(alignment.del, deletion), unlist)
          alignment.del.valid <- sapply(alignment.del, function(x) {length(x) > 0})
          alignment.ins.valid <- sapply(alignment.ins, function(x) {length(x) > 0})
          full.replace <- (alignment.del.valid & alignment.ins.valid & alt.width == ref.width)
          alignment.del.valid <- alignment.del.valid & !full.replace
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
                           with.indels = FALSE, max.mismatch = 0)
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
  return(list(fsnplist = fsnplist,
              ref.seq = snp.sequence.ref,
              alt.seq = snp.sequence.alt))
}

## An evaluator function for SNP effect

 varEff <- function(allelR, allelA) {
  score <- allelA - allelR
  if (abs(score) >= 0.7) {
    return(list(score = score, effect = "strong"))
  } else if (abs(score) > 0.4) {
    return(list(score = score, effect = "weak"))
  } else {
    return(list(score = score, effect = "neut"))
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


scoreSeqWindows <- function(ppm, seq) {
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
  scores_rc <- t(reverseComplementMotif(ppm)[seq, ])[ranges]
  scores <- split(scores, ceiling(seq_along(scores)/ppm.width))
  scores_rc <- split(scores_rc, ceiling(seq_along(scores_rc)/ppm.width))
  res <- vapply(Map(function(x, y) {matrix(data = c(x, y), nrow = 2,
                                           byrow = TRUE, dimnames = list(c(1, 2)))},
                    scores, scores_rc),
                rowSums,
                numeric(2))
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

  snp.sequence.alt <- fsnplist$alt.seq
  snp.sequence.ref <- fsnplist$ref.seq
  fsnplist <- fsnplist$fsnplist

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
    res.el$scoreRef <- as.numeric(NA)
    res.el$scoreAlt <- as.numeric(NA)
    if (filterp) {
      res.el$Refpvalue <- as.numeric(NA)
      res.el$Altpvalue <- as.numeric(NA)
    }
    if (ref.len > 1 | alt.len > 1) {
      res.el$altPos <- as.numeric(NA)
      res.el$alleleDiff <- as.numeric(NA)
      res.el$alleleEffectSize <- as.numeric(NA)
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
      seq.start <- min(alt.loc)
      seq.len <- length(alt.loc)
      alt.range <- ref.range <- (k - (ncol(pwm) - seq.start)):(k + ncol(pwm) + seq.start + seq.len - 2)
      if (!show.neutral & identical(snp.ref[ref.range], snp.alt[alt.range])) next()
      seq.remove <- ref.len - alt.len
      if (seq.remove < 0) {
        ref.range <- ref.range[1:(length(ref.range) + seq.remove)]
      } else {
        alt.range <- alt.range[1:(length(alt.range) - seq.remove)]
      }
      ref.windows <- scoreSeqWindows(ppm = pwm, seq = snp.ref[ref.range])
      alt.windows <- scoreSeqWindows(ppm = pwm, seq = snp.alt[alt.range])
      pass.effect <- ifelse(filterp,
                            any(alt.windows > thresh) | any(ref.windows > thresh),
                            any(((alt.windows - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1]) > thresh)) |
                                any((ref.windows - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1]) > thresh))
      if (pass.effect) {
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
          allelR <- ref.windows[hit.ref$strand, hit.ref$window]
          allelA <- alt.windows[hit.alt$strand, hit.alt$window]
          scorediff <- varEff(allelR, allelA)
          effect <- scorediff$effect
          score <- scorediff$score
          #ref.pos <- ref.range[(ncol(pwm) - (seq.start - 1)):((ncol(pwm) + ref.len - seq.start))]
          ref.pos <- k:(k + nchar(result$REF) - 1L)
          #alt.pos <- alt.range[(ncol(pwm) - (seq.start - 1)):((ncol(pwm) + alt.len - seq.start))]
          alt.pos <- k:(k + nchar(result$ALT) - 1L)
          if ((effect == "neut" & show.neutral) | effect != "neut") {
            res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                         snp.ref, snp.alt,
                                                         ref.pos, alt.pos,
                                                         hit.ref, hit.alt,
                                                         ref.windows, alt.windows,
                                                         score, effect, len,
                                                         k, pwm, calcp = filterp)
          }
        } else {
          snp.pos <- k:(k + nchar(result$REF) - 1L)
          allelR <- ref.windows[hit.ref$strand, hit.ref$window]
          allelA <- alt.windows[hit.alt$strand, hit.alt$window]
          scorediff <- varEff(allelR, allelA)
          effect <- scorediff$effect
          score <- scorediff$score
          if ((effect == "neut" & show.neutral) | effect != "neut") {
            res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                         snp.ref, snp.alt,
                                                         snp.pos, snp.pos,
                                                         hit.ref, hit.alt,
                                                         ref.windows, alt.windows,
                                                         score, effect, len,
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

#' @importFrom matrixStats colRanges
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
  alt_loc <- range(mresult$ALT_loc)
  ref_start <- (1 - alt_loc[[1]])
  ref_start <- ifelse(ref_start <= 0, ref_start - 1, ref_start)
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
  pwmrange <- colSums(colRanges(pwm[-5,]))
  mresult[["scoreRef"]] <- ref.windows[hit.ref$strand, hit.ref$window]
  mresult[["scoreAlt"]] <- alt.windows[hit.alt$strand, hit.alt$window]
  mresult[["pctRef"]] <- (mresult[["scoreRef"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
  mresult[["pctAlt"]] <- (mresult[["scoreAlt"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
  if (calcp) {
    mresult[["Refpvalue"]] <- NA
    mresult[["Altpvalue"]] <- NA
  }
  mresult[["alleleDiff"]] <- score
  mresult[["effect"]] <- effect
  mresult[["alleleEffectSize"]] <- score/pwmrange[[2]]
  mcols(result) <- mresult
  return(result)
}

#' @importFrom matrixStats colMaxs colMins
preparePWM <- function(pwmList,
                       filterp,
                       bkg,
                       scoreThresh,
                       method = "default") {

  bkg <- bkg[c('A', 'C', 'G', 'T')]

  scounts <- as.integer(mcols(pwmList)$sequenceCount)
  scounts[is.na(scounts)] <- 20L
  pwmList.pc <- Map(function(pwm, scount) {
    pwm <- (pwm * scount + 0.25)/(scount + 1)
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
      pwm <- log(pwm) - log(b)
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
#'  \item{REF}{the reference allele for the variant}
#'  \item{ALT}{the alternate allele for the variant}
#'  \item{snpPos}{the coordinates of the variant}
#'  \item{motifPos}{The position of the motif relative the the variant}
#'  \item{geneSymbol}{the geneSymbol corresponding to the TF of the TF binding motif}
#'  \item{dataSource}{the source of the TF binding motif}
#'  \item{providerName, providerId}{the name and id provided by the source}
#'  \item{seqMatch}{the sequence on the 5' -> 3' direction of the "+" strand
#'  that corresponds to DNA at the position that the TF binding motif was found.}
#'  \item{pctRef}{The score as determined by the scoring method, when the sequence contains the reference variant allele, normalized to a scale from 0 - 1. If \code{filterp = FALSE},
#'  this is the value that is thresholded.}
#'  \item{pctAlt}{The score as determined by the scoring method, when the sequence contains the alternate variant allele, normalized to a scale from 0 - 1. If \code{filterp = FALSE},
#'  this is the value that is thresholded.}
#'  \item{scoreRef}{The score as determined by the scoring method, when the sequence contains the reference variant allele}
#'  \item{scoreAlt}{The score as determined by the scoring method, when the sequence contains the alternate variant allele}
#'  \item{Refpvalue}{p-value for the match for the pctRef score, initially set to \code{NA}. see \code{\link{calculatePvalue}} for more information}
#'  \item{Altpvalue}{p-value for the match for the pctAlt score, initially set to \code{NA}. see \code{\link{calculatePvalue}} for more information}
#'  \item{alleleRef}{The proportional frequency of the reference allele at position \code{motifPos} in the motif}
#'  \item{alleleAlt}{The proportional frequency of the alternate allele at position \code{motifPos} in the motif}
#'  \item{altPos}{the position, relative to the reference allele, of the alternate allele}
#'  \item{alleleDiff}{The difference between the score on the reference allele and the score on the alternate allele}
#'  \item{alleleEffectSize}{The ratio of the \code{alleleDiff} and the maximal score of a sequence under the PWM}
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
#' @import parallel
#' @importFrom parallel clusterEvalQ
#' @importFrom BiocParallel bplapply
#' @importFrom stringr str_length str_trim
#' @export
motifbreakR <- function(snpList, pwmList, threshold = 0.85, filterp = FALSE,
                        method = "default", show.neutral = FALSE, verbose = FALSE,
                        bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                        BPPARAM = bpparam()) {
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
  if (!(is(BPPARAM, "MulticoreParam")|is(BPPARAM, "SerialParam"))) {
    bpstart(BPPARAM)
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  ##
  ## Genome Setup
  genome.package <- attributes(snpList)$genome.package
  if (requireNamespace(eval(genome.package), quietly = TRUE, character.only = TRUE)) {
    genome.bsgenome <- eval(parse(text = paste(genome.package, genome.package, sep = "::")))
  } else {
    stop(paste0(eval(genome.package), " is the genome selected for this snp list and \n",
                "is not present on your environment. Please load it and try again."))
  }
  ##

  pwms <- preparePWM(pwmList = pwmList, filterp = filterp,
                     scoreThresh = threshold, bkg = bkg,
                     method = method)

  k <- max(sapply(pwms$pwmList, ncol))

  snpList <- prepareVariants(fsnplist = snpList,
                             genome.bsgenome = genome.bsgenome,
                             max.pwm.width = k)

  snpList_cores <- split(as.list(rep(names(snpList), times = cores)), 1:cores)
  for (splitr in seq_along(snpList)) {
    splitcores <- sapply(suppressWarnings(split(snpList[[splitr]], 1:cores)), list)
    for (splitcore in seq_along(snpList_cores)) {
      snpList_cores[[splitcore]][[splitr]] <- splitcores[[splitcore]]
      names(snpList_cores[[splitcore]])[splitr] <- names(snpList)[splitr]
    }
  }
  snpList <- snpList_cores; rm(snpList_cores)

  x <- bplapply(snpList, scoreSnpList,
                pwmList = pwms$pwmList, threshold = pwms$pwmThreshold,
                pwmList.pc = pwms$pwmListPseudoCount, pwmRanges = pwms$pwmRange,
                method = method, bkg = bkg, show.neutral = show.neutral,
                verbose = ifelse(cores == 1, verbose, FALSE), genome.bsgenome = genome.bsgenome,
                filterp = filterp, BPPARAM = BPPARAM)

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
    snpList <- unlist(GRangesList(lapply(snpList, `[[`, "fsnplist")), use.names = FALSE)
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
#' @param background Numeric Vector; the background probabilities of the nucleotides
#' @param granularity Numeric Vector; the granularity to which to round the PWM,
#'  larger values compromise full accuracy for speed of calculation. A value of
#'  \code{NULL} does no rounding.
#' @param BPPARAM a BiocParallel object see \code{\link[BiocParallel]{register}}
#'   and see \code{getClass("BiocParallelParam")} for additional parameter
#'   classes.  Try \code{BiocParallel::registered()} to see what's available and
#'   for example \code{BiocParallel::bpparam("SerialParam")} would allow serial
#'   evaluation.
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
#' rs1006140 <- example.results[example.results$SNP_id %in% "rs1006140"]
#' # low granularity for speed; 1e-6 or 1e-7 recommended for accuracy
#' rs1006140 <- calculatePvalue(rs1006140, BPPARAM=BiocParallel::SerialParam(), granularity = 1e-4)
#'
# #' @importFrom qvalue qvalue
#'
#' @export
calculatePvalue <- function(results,
                            background = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                            granularity = NULL,
                            BPPARAM = BiocParallel::SerialParam()) {

  ## Cluster / MC setup
  if (.Platform$OS.type == "windows" && inherits(BPPARAM, "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n",
                   "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpnworkers(BPPARAM)
  num.res <- length(results)
  if (num.res < cores) {
    cores <- num.res
  }
  if (!(is(BPPARAM, "MulticoreParam")|is(BPPARAM, "SerialParam"))) {
    bpstart(BPPARAM)
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  if(!("scoreRef" %in% names(mcols(results)))) {
    stop('incorrect results format; please rerun analysis with filterp=TRUE')
  } else {
    pwmListmeta <- mcols(attributes(results)$motifs, use.names=TRUE)
    pwmList <- attributes(results)$scoremotifs
    if(!is.null(granularity)) {
      pwmList <- lapply(pwmList, function(x, g) {x <- floor(x/g)*g; return(x)}, g = granularity)
    }
    results_sp <- split(results, 1:length(results))
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
    ## Cluster / MC cleanup
    if (inherits(pvalues, "try-error")) {
      if (is(BPPARAM, "SnowParam")) {
        bpstop(BPPARAM)
      }
      stop(attributes(pvalues)$condition)
    }
    pvalues.df <- base::do.call("rbind", c(pvalues, make.row.names = FALSE))
    results$Refpvalue <- pvalues.df[, "ref"]
    results$Altpvalue <- pvalues.df[, "alt"]

    if (is(BPPARAM, "SnowParam")) {
      bpstop(BPPARAM)
    }

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
    #if (class(.ele) != "pfm")
    if (!is(.ele, 'pfm'))
      stop("pfms must be a list of class pfm")
  })
  assign("tmp_motifStack_symbolsCache", list(), pos = ".GlobalEnv")
  # grid.newpage()
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
#' @param altAllele Character; The default value of \code{NULL} uses the first (or only)
#'   alternative allele for the SNP to be plotted.
#' @seealso See \code{\link{motifbreakR}} for the function that produces output to be
#'   visualized here, also \code{\link{snps.from.rsid}} and \code{\link{snps.from.file}}
#'   for information about how to generate the input to \code{\link{motifbreakR}}
#'   function.
#' @details \code{plotMB} produces output showing the location of the SNP on the
#'   chromosome, the surrounding sequence of the + strand, the footprint of any
#'   motif that is disrupted by the SNP or SNV, and the DNA sequence motif(s).
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
#' @importFrom Gviz IdeogramTrack SequenceTrack GenomeAxisTrack HighlightTrack
#'   AnnotationTrack plotTracks
#' @export
plotMB <- function(results, rsid, reverseMotif = TRUE, effect = c("strong", "weak"), altAllele = NULL) {
  motif.starts <- sapply(results$motifPos, `[`, 1)
  motif.starts <- start(results) + motif.starts
  motif.starts <- order(motif.starts)
  results <- results[motif.starts]
  g <- genome(results)[[1]]
  result <- results[results$SNP_id %in% rsid]
  if(is.null(altAllele)) {
    altAllele <- result$ALT[[1]]
  }
  result <- result[result$ALT == altAllele]
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
  results_motifs <- paste0(result$providerId, result$providerName)
  list_motifs <- paste0(mcols(pwmList)$providerId, mcols(pwmList)$providerName)
  pwms <- pwmList <- pwmList[match(results_motifs, list_motifs)]
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
  } else {
    for (pwm.i in seq_along(pwms)) {
      pwm.name <- names(pwms[pwm.i])
      pwm.id <- mcols(pwms[pwm.name, ])$providerId
      pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
      doRev <- as.logical(strand(result[result$providerId == pwm.id & result$providerName == pwm.name.f, ]) == "-")
      if (doRev) {
        pwm <- pwms[[pwm.i]]
        pwm <- pwm[, rev(1:ncol(pwm))]
        pwms[[pwm.i]] <- pwm
        names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], "-:r")
      }
    }
  }
  pwms <- lapply(names(pwms), function(x, pwms=pwms) {new("pfm", mat = pwms[[x]],
                                                          name = x)}, pwms)
  pwms <- DNAmotifAlignment.2snp(pwms, result)
  pwmwide <- max(sapply(pwms, function(x) { ncol(x@mat)}))

  markerStart <- result$motifPos[[1]][1]
  if (markerStart > 0) {
    markerEnd <- length(result$altPos[[1]]) + 1
    markerEnd <- markerEnd - markerStart
    markerStart <- 1
  } else {
    markerStart <- -1 * markerStart
    markerEnd <- markerStart + length(result$altPos[[1]])
    if (result$varType[[1]] %in% c("Other", "SNV")) {
      markerStart <- markerStart + 1
    }
  }
  varType <- result$varType[[1]]
  varType <- switch(varType,
                    Deletion = "firebrick",
                    Insertion = "springgreen4",
                    Other = "gray13")
  markerRect <- new("marker", type = "rect",
                    start = markerStart,
                    stop = markerEnd,
                    gp = gpar(lty = 2,
                              fill = NA,
                              lwd = 3,
                              col = varType))
  for (pwm.i in seq_along(pwms)) {
    pwms[[pwm.i]]@markers <- list(markerRect)
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

  #altseq <- replaceLetterAt(altseq, at = wherereplace, letter = rep.int("N", sum(wherereplace)))
  #altseq <- replaceLetterAt(altseq, at = !wherereplace, letter = result$ALT[[1]])
  histart <- start(result[1]) + min(result[1]$altPos[[1]]) - 2
  histart <- ifelse(result[1]$varType %in% c("Other", "SNV"), histart + 1, histart)
  hiend <- start(result[1]) + min(result[1]$altPos[[1]]) - 2 + length(result[1]$altPos[[1]])
  hiT <- HighlightTrack(trackList = list(seqT, seqAltT),
                        start = histart,
                        end = hiend,
                        chromosome = chromosome)

  selectingfun <- selcor
  detailfun <- addPWM.stack

  motif_ids <- names(pwmList)
  names(motif_ids) <- mcols(pwmList)$providerName

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
  presult <- result
  strand(presult) <- "*"
  pres_cols <- DataFrame(feature = ifelse(!is.na(result$motifPos),
                                          paste(result$geneSymbol, "motif", sep = "_"), ""),
                         group = result$providerName,
                         id = motif_ids)
  presult <- GRanges(seqnames = seqnames(result[1]),
                     ranges = ranges(result))
  mcols(presult) <- pres_cols

  motifT <- AnnotationTrack(presult,
                            fun = detailfun,
                            detailsFunArgs = list(pwm_stack = pwms),
                            name = names(result)[[1]],
                            selectFun = selectingfun,
                            reverseStacking = FALSE,
                            stacking = "squish")

  if (exists("axisT")) {
    track_list <- list(ideoT, motifT, hiT, axisT)
  } else {
    track_list <- list(ideoT, motifT, hiT)
  }
  plotTracks(track_list, from = from, to = to, showBandId = TRUE,
             cex.main = 0.8, col.main = "darkgrey",
             add53 = TRUE, labelpos = "below", chromosome = chromosome, #groupAnnotation = "id",
             fontcolor.item="black",
             collapse = FALSE, min.width = 1, featureAnnotation = "feature", cex.feature = 0.8,
             details.size = 0.85, detailsConnector.pch = NA, detailsConnector.lty = 0,
             shape = "box", cex.group = 0.8, fonts = c("sans", "Helvetica"))
  return(invisible(NULL))
}

