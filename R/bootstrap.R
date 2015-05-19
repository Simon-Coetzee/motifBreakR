#' @importFrom Biostrings DNAString
bootstrap.for.empirical.pvalue <- function(pwmList, method = c("default", "ic", "log"), BPPARAM) {
  seq <- sample(c("A", "C", "G", "T"), 11000, replace = TRUE)
  seq <- DNAString(paste(seq, collapse = ""))
  seq.rc <- reverseComplement(seq)
  seq <- strsplit(as.character(seq), "")[[1]]
  seq.rc <- strsplit(as.character(seq.rc), "")[[1]]
  results <- list()
  bkg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
  for(meth in method) {
    scounts <- as.integer(mcols(pwmList)$sequenceCount)
    scounts[is.na(scounts)] <- 20L
    pwmList.pc <- Map(function(pwm, scount) {
                        pwm <- (pwm * scount + 0.25)/(scount + 1)
                      }, pwmList, scounts)
    if(meth == "ic") {
      pwmOmegas <- lapply(pwmList.pc, function(pwm, b=bkg) {
                            omegaic <- colSums(pwm * log2(pwm/b))
                          })
    }
    if(meth == "default") {
      pwmOmegas <- lapply(pwmList.pc, function(pwm) {
                            omegadefault <- colMaxs(pwm) - colMins(pwm)
                          })
    }
    if(meth == "log") {
      pwmOmegas <- 1
    }
    pwmRanges <- Map(function(pwm, omega) {
                       x <- colSums(colRanges(pwm) * omega)
                       return(x)
                     }, pwmList.pc, pwmOmegas)
    motifin <- Map(list, pwmList.pc, pwmOmegas, pwmRanges, names(pwmList))
    result <- bplapply(motifin, function(ppm, me, ...) {
                       b <- c(A=0.25, C=0.25, G=0.25, T=0.25)
                       pwm <- ppm[[1]]
                       scorer <- scoreAllWindows(snp.seq = seq,
                                                 snp.seq.rc = seq.rc,
                                                 pwm = pwm,
                                                 method = me,
                                                 from = 1,
                                                 to = length(seq) - 1 - ncol(pwm),
                                                 bkg = b,
                                                 pwm.range = ppm[[3]],
                                                 omega = ppm[[2]],
                                                 pvalue = TRUE)
                       scorer <- as.numeric(scorer[1,])
                       return(scorer)
                     }, me = meth, seq, seq.rc, BPPARAM=BPPARAM)
    result <- lapply(result, sample, 10000)
    result <- do.call(cbind, result)
    result <- list(result)
    names(result) <- meth
    results <- c(results, result)
    rm(result)
    gc()
  }
  return(results)
}
