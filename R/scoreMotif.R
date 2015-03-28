## scoreMotif.R this module will score an arbitrary list of SNPs for motif
## disruption using a buffet of menu options for some sets of transcription
## factors

## optional: produce random SNP sets of equivalent size to the index SNP list to
## simulate background for enrichment calculations.

require("BSgenome.Hsapiens.UCSC.hg19") ## DNAString, DNAStringSet, PDict, matchPDict
require("MotifDb") ## dataset

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
  f_b <- ppm
  i <- ncol(f_b)
  p_b <- matrix(rep(bkg, i), nrow = 4)
  row.names(p_b) <- row.names(f_b)
  IC_i <- f_b * t(apply(f_b/p_b + 1e-07, 1, log2))
  IC <- apply(IC_i, 2, sum)
  IC
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
  if (method == "log") {
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
    ic_pwm <- defaultIC(pwm, bkg)
    vscores <- PWMlimit * ic_pwm
  }
  sum(vscores)
}


pctScore <- function(score, min.score, max.score) {
  (score - min.score)/(max.score - min.score)
}

## motif scorer, using MotifDb objects returns vector of probabilities based on
## the PPM
#' importFrom compiler compfun
scoreMotif <- function(snp.seq, ppm, len, method = "default",
                       ag = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                       offset = 1) {
  snp.seq <- snp.seq[offset:(offset + len - 1)]
  ## diag code
  position.probs <- c(ppm[snp.seq, ])[1 + 0L:(len - 1L) * (len + 1)]
  if (method == "default") {
    return(position.probs)
  } else {
    position.probs <- position.probs/ag[snp.seq]
    return(position.probs)
  }
}
scoreMotif2 <- cmpfun(scoreMotif)
require(compiler)

## create an all purpose scoring algorithm 'w' for weighted
wScore <- function(snp.seq, ppm, pwm.omega, offset, method = "default", bkg = NULL, 
                   len) {
  if (is.null(bkg)) {
    bkg <- c(0.25, 0.25, 0.25, 0.25)
  } else {
    bkg <- bkg[c(1, 2, 2, 1)]
  }
  names(bkg) <- c("A", "C", "G", "T")
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
    ic_ppm <- defaultIC(ppm, bkg)
    returnScore <- sum(vscores * ic_ppm)
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

scoreSnpList <- function(fsnplist, pwmList, method = "default", bkg = NULL, threshold = 0.8, 
                         show.neutral = FALSE, verbose = FALSE) {
  if (!(Reduce("&", width(fsnplist$REF) == 1)) || !(Reduce("&", width(fsnplist$ALT) == 
                                                             1))) {
    drop.snp.i <- !(width(fsnplist$REF) == 1) & (width(fsnplist$ALT) == 1)
    drop.snp <- names(fsnplist[!drop.snp.i])
    fsnplist <- fsnplist[drop.snp.i]
    warning(paste(drop.snp, "not valid for motif analysis...skipping", sep = " "))
  }
  genome.package <- attributes(fsnplist)$genome.package
  genome.bsgenome <- eval(parse(text = genome.package))
  require(eval(genome.package), character.only = TRUE)
  k <- max(sapply(pwmList, ncol))
  ## Get Reference Sequence with k bp flanking
  snp.sequence.ref <- getSeq(genome.bsgenome, promoters(fsnplist, upstream = k - 
                                                          1, downstream = k))
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
    return(rep.int(x, l))
  })
  pwmList.omega <- lapply(pwmList, defaultOmega)
  pwmList.max <- lapply(pwmList, limitPwm, limitFun = max, method = method, bkg = bkg)
  pwmList.min <- lapply(pwmList, limitPwm, limitFun = min, method = method, bkg = bkg)
  pwmList.bot <- mapply("-", pwmList.max, pwmList.min)
  strand.opt <- c("+", "-")
  set.sig <- rep(NA, length(resultSet))
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
    pwm.sig <- rep(NA, length(res.el))
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
        pwm.sig[[pwm.i]] <- FALSE
      } else {
        if (hit.ref[["strand"]] != 0L) {
          hit <- hit.ref
          snp.pos <- len - hit[["window"]] + 1
          result <- res.el[pwm.i]
          if (hit[["strand"]] == 1L) {
            allelR <- pwm[as.character(result$REF), snp.pos]
            allelA <- pwm[as.character(result$ALT), snp.pos]
          } else {
            allelR <- pwm[as.character(reverse(result$REF)), snp.pos]
            allelA <- pwm[as.character(reverse(result$ALT)), snp.pos]
          }
          effect <- snpEff(allelR, allelA)
          if (effect == "neut") {
            if (show.neutral) {
              res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos, 
                                             hit, ref.windows, alt.windows, allelR, allelA, effect, len, 
                                             k, pwm)
              pwm.sig[[pwm.i]] <- TRUE
            } else {
              pwm.sig[[pwm.i]] <- FALSE
            }
          } else {
            res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos, 
                                           hit, ref.windows, alt.windows, allelR, allelA, effect, len, 
                                           k, pwm)
            pwm.sig[[pwm.i]] <- TRUE
          }
        }
        if (hit.alt[["strand"]] != 0L && hit.alt[["window"]] != hit.ref[["window"]]) {
          hit <- hit.alt
          snp.pos <- len - hit[["window"]] + 1
          result <- res.el[pwm.i]
          if (hit[["strand"]] == 1L) {
            allelR <- pwm[as.character(result$REF), snp.pos]
            allelA <- pwm[as.character(result$ALT), snp.pos]
          } else {
            allelR <- pwm[as.character(reverse(result$REF)), snp.pos]
            allelA <- pwm[as.character(reverse(result$ALT)), snp.pos]
          }
          effect <- snpEff(allelR, allelA)
          if (effect == "neut") {
            if (show.neutral) {
              res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos, 
                                             hit, ref.windows, alt.windows, allelR, allelA, effect, len, 
                                             k, pwm)
              pwm.sig[[pwm.i]] <- TRUE
            } else {
              pwm.sig[[pwm.i]] <- FALSE
            }
          } else {
            res.el[pwm.i] <- updateResults(res.el[pwm.i], snp.ref, snp.pos, 
                                           hit, ref.windows, alt.windows, allelR, allelA, effect, len, 
                                           k, pwm)
            pwm.sig[[pwm.i]] <- TRUE
          }
        }
      }
    }
    if (verbose) {
      message(names(res.el[1]))
      message(paste0("number of broken motifs = ", sum(pwm.sig)))
    }
    res.el <- res.el[pwm.sig]
    if (length(res.el) > 0) {
      resultSet[[snp.map.i]] <- res.el
      set.sig[[snp.map.i]] <- TRUE
    } else {
      set.sig[[snp.map.i]] <- FALSE
    }
  }
  resultSet <- resultSet[set.sig]
  if (length(resultSet) < 1) {
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist), 
                    "with 0 potentially disruptive matches to", length(unique(resultSet$geneSymbol)), 
                    "of", length(pwmList), "motifs."))
    }
    return(NULL)
  } else {
    resultSet <- unlist(GRangesList(resultSet), use.names = F)
    mcols(resultSet) <- mcols(resultSet)[, -(1:2)]
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist), 
                    "with", length(resultSet), "potentially disruptive matches to", length(unique(resultSet$geneSymbol)), 
                    "of", length(pwmList), "motifs."))
    }
    return(resultSet)
  }
}


updateResults <- function(result, snp.seq, snp.pos, hit, ref.windows, alt.windows, 
                          allelR, allelA, effect, len, k, pwm) {
  ## if(snp.pos == 6) stop('here')
  strand.opt <- c("+", "-")
  strand <- (-hit[["strand"]] + 3) / 2
  strand(result) <- strand.opt[strand]
  result$snpPos <- start(result)
  result$motifPos <- snp.pos
  if (strand == 2) {
    matchs <- snp.seq[(k - hit[["window"]] + 1):(k - hit[["window"]] + len)]
    matchs[-hit[["window"]]] <- tolower(matchs[-hit[["window"]]])
    matchs <- paste(matchs, collapse = "")
    result$seqMatch <- str_pad(matchs, width = k + hit[["window"]] - 1, side = "right")
    ## pwm pos instead of snp pos
    start(result) <- start(result) - hit[["window"]] + 1
    end(result) <- end(result) - hit[["window"]] + len
  } else {
    matchs <- snp.seq[(k - len + hit[["window"]]):(k + hit[["window"]] - 1)]
    matchs[-snp.pos] <- tolower(matchs[-snp.pos])
    matchs <- paste(matchs, collapse = "")
    result$seqMatch <- str_pad(matchs, width = k + snp.pos - 1, side = "right")
    start(result) <- start(result) - snp.pos + 1
    end(result) <- end(result) - snp.pos + len
  }
  result$pctRef <- ref.windows[strand, hit[["window"]]]
  result$pctAlt <- alt.windows[strand, hit[["window"]]]
  result$alleleRef <- allelR
  result$alleleAlt <- allelA
  result$effect <- effect
  return(result)
}

motifbreakR <- function(snpList, pwmList, threshold, method = "default", bkg = NULL, 
                        show.neutral = FALSE, verbose = FALSE, cores = 1L) {
  require(parallel)
  num.snps <- length(snpList)
  genome.package <- attributes(snpList)$genome.package
  snpList <- sapply(suppressWarnings(split(snpList, 1:cores)), list)
  x <- mclapply(snpList, scoreSnpList, pwmList = pwmList, threshold = threshold, 
                method = method, bkg = bkg, show.neutral = show.neutral, verbose = ifelse(cores == 
                                                                                            1, verbose, FALSE), mc.cores = cores, mc.preschedule = FALSE)
  drops <- sapply(x, is.null)
  x <- x[!drops]
  if (length(x) > 1) {
    x <- unlist(GRangesList(unname(x)))
    x <- x[order(match(names(x), names(snpList)), x$geneSymbol), ]
    attributes(x)$genome.package <- genome.package
    attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% unique(x$providerId), 
                                    ]
  } else {
    if (length(x) == 1L) {
      x <- x[[1]]
      attributes(x)$genome.package <- genome.package
      attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% unique(x$providerId), 
                                      ]
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


addPWM <- function(identifier, ...) {
  motif.name <- identifier
  ps.file <- paste(tempdir(), motif.name, sep = "/")
  PostScriptTrace(ps.file, paste0(ps.file, ".xml"))
  motif.figure <- readPicture(paste0(ps.file, ".xml"))
  grid.picture(motif.figure[-1], distort = FALSE)
}

MBplot <- function(results, rsid, reverseMotif = TRUE) {
  g <- genome(results)[[1]]
  result <- results[names(results) %in% rsid]
  chromosome <- as.character(seqnames(result))[[1]]
  genome.package <- attributes(result)$genome.package
  genome.bsgenome <- eval(parse(text = genome.package))
  from <- min(start(result)) - 3
  to <- max(end(result)) + 4
  temp.dir <- tempdir()
  pwmList <- attributes(result)$motifs
  for (pwm.id in names(pwmList)) {
    postscript(paste(temp.dir, pwm.id, sep = "/"), width = 10, height = 3, horizontal = FALSE, 
               fonts = c("sans", "Helvetica"))
    pwm <- pwmList[[pwm.id]]
    if (Reduce("|", result$providerId %in% mcols(pwmList[pwm.id])$providerId)) {
      if (reverseMotif && (as.character(strand(result)[result$providerId %in% 
                                                       mcols(pwmList[pwm.id])$providerId]) == "-")) {
        pwm <- pwm[, rev(1:ncol(pwm))]
        rownames(pwm) <- c("T", "G", "C", "A")
      }
    }
    p <- new("pfm", mat = pwm, name = mcols(pwmList[pwm.id])$providerId)
    plot(p)
    dev.off()
  }
  ideoT <- IdeogramTrack(genome = g, chromosome = chromosome)
  seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", "auto"))
  axisT <- GenomeAxisTrack(exponent = 0)
  hiT <- HighlightTrack(seqT, start = result$snpPos[[1]], end = result$snpPos[[1]], 
                        chromosome = chromosome)
  if (Reduce("|", result$providerId %in% mcols(pwmList[pwm.id])$providerId)) {
    if (reverseMotif && (as.character(strand(result)[result$providerId %in% mcols(pwmList[pwm.id])$providerId]) == 
                           "-")) {
      motifT <- AnnotationTrack(result, id = names(pwmList)[mcols(pwmList)$providerId %in% 
                                          result$providerId], fun = addPWM, group = result$providerId, feature = paste0("snp@", 
                                                                                                         result$motifPos), name = names(result)[[1]])
    }
  } else {
    motifT <- AnnotationTrack(result, id = names(pwmList)[mcols(pwmList)$providerId %in% 
                                        result$providerId], fun = addPWM, group = result$providerId, feature = paste0("snp@", 
                                                                                                       result$motifPos), name = names(result)[[1]])
  }
  plotTracks(list(ideoT, motifT, hiT, axisT), from = from, to = to, showBandId = TRUE, 
             add53 = TRUE, labelpos = "below", chromosome = chromosome, groupAnnotation = "group", 
             collapse = FALSE, min.width = 1, featureAnnotation = "feature", cex.feature = 0.8, 
             details.size = 0.5, detailsConnector.pch = NA, shape = "box", cex.group = 0.8, 
             fonts = c("sans", "Helvetica"))
}



