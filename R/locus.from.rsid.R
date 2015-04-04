#' @importFrom BSgenome snpid2grange snplocs
#' @importFrom Biostrings DNAStringSet
#' @export
snps.from.rsid <- function(rsid = NULL, dbSNP = NULL,
                           search.genome = NULL) {
  if (is.null(rsid)) {
    stop("no RefSNP IDs have been found, please include RefSNP ID numbers")
  }
  if (class(dbSNP) != "SNPlocs") {
    stop(paste0("dbSNP argument was not provided with a valid SNPlocs object.\n",
                "Please run availible.SNPs() to check for availble SNPlocs"))
  }
  if (class(search.genome) != "BSgenome") {
    stop(paste0("search.genome argument was not provided with a valid BSgenome object.\n",
                "Run availible.genomes() and choose the appropriate BSgenome object"))
  }
  if (Reduce("&", !grepl("rs", rsid))) {
    bad.names <- rsid[!grepl("rs", rsid)]
    stop(paste(paste(bad.names, collapse = " "), "are not rsids, perhaps you want to import your snps from a bed file with snps.from.bed()?"))
  }
  rsid.grange <- snpid2grange(dbSNP, rsid)
  rsid.grange <- change.to.search.genome(rsid.grange, search.genome)
  rsid.refseq <- getSeq(search.genome, rsid.grange)
  rsid.grange$UCSC.reference <- as.character(rsid.refseq)
  rsid.grange <- sapply(split(rsid.grange, rsid.grange$RefSNP_id), function(snp) {
    alt.allele <- determine.allele.from.ambiguous(snp$alleles_as_ambig, snp$UCSC.reference)
    if (length(alt.allele) > 1L) {
      snp <- do.call("c", replicate(length(alt.allele), snp))
      snp$UCSC.alternate <- alt.allele
      names(snp) <- paste(gsub("^", "rs", snp$RefSNP_id), alt.allele, sep = ":")
    } else {
      snp$UCSC.alternate <- alt.allele
      names(snp) <- gsub("^", "rs", snp$RefSNP_id)
    }
    return(snp)
  })
  rsid.grange <- unlist(do.call("GRangesList", rsid.grange), use.names = FALSE)
  colnames(mcols(rsid.grange)) <- c("RefSNP_id", "alleles_as_ambig", "REF", "ALT")
  rsid.grange$REF <- DNAStringSet(rsid.grange$REF)
  rsid.grange$ALT <- DNAStringSet(rsid.grange$ALT)
  rsid.grange$alleles_as_ambig <- DNAStringSet(rsid.grange$alleles_as_ambig)
  colnames(mcols(rsid.grange))[1] <- "SNP_id"
  attributes(rsid.grange)$genome.package <- attributes(search.genome)$pkgname
  return(rsid.grange)
}

determine.allele.from.ambiguous <- function(ambiguous.allele, known.allele) {
  neucleotide.ambiguity.code <- list(Y = c("C", "T"), R = c("A", "G"), W = c("A", "T"),
                                     S = c("G", "C"), K = c("T", "G"), M = c("C", "A"),
                                     D = c("A", "G", "T"), V = c("A", "C", "G"),
                                     H = c("A", "C", "T"), B = c("C", "G", "T"),
                                     N = c("A", "C", "G", "T"))
  specnac <- neucleotide.ambiguity.code[[ambiguous.allele]]
  unknown.allele <- specnac[-grep(known.allele, specnac)]
  return(unknown.allele)
}


#change.to.ucsc.genome <- function(granges.object, search.genome = BSgenome.Hsapiens.UCSC.hg19) {
#  if (Reduce("&", !is.na(genome(granges.object)))) {
#    if (Reduce("&", genome(granges.object) == "hg19")) {
#      return(granges.object)
#    }
#  }
#  if (Reduce("&", is.na(genome(granges.object))) && length(seqinfo(granges.object)) !=
#      25) {
#    xome.value <- str_extract(seqlevels(granges.object), "[0-9]|1[0-9]|2[0-9]|3[0-9]|4[0-9]|5[0-9]|6[0-9]|7[0-9]|8[0-9]|9[0-9]|M")
#    positions <- sapply(paste0("chr", xome.value, "$"), grep, seqnames(seqinfo(search.genome)))
#    new2oldmap <- rep(NA, length(seqinfo(search.genome)))
#    new2oldmap[positions] <- 1:length(positions)
#  } else {
#    new2oldmap <- c(1:25, rep(NA, length(seqinfo(search.genome)) - length(seqinfo(granges.object))))
#  }
#  seqinfo(granges.object, new2old = new2oldmap) <- seqinfo(search.genome)
#  return(granges.object)
#}

#' @importFrom stringr str_extract
#' @importFrom GenomeInfoDb seqinfo genome
#' @importFrom GenomeInfoDb seqnames seqlevels seqinfo<-
change.to.search.genome <- function(granges.object, search.genome) {
  if (Reduce("&", !is.na(genome(granges.object)))) {
    if (identical(genome(granges.object), genome(search.genome))) {
      return(granges.object)
    }
  }
  if(isTRUE(all.equal(seqlevels(granges.object), seqlevels(search.genome)))) {
    seqinfo(granges.object) <- seqinfo(search.genome)
  } else {
    if(seqlevelsStyle(granges.object) != seqlevelsStyle(search.genome)) {
      seqlevelsStyle(granges.object) <- seqlevelsStyle(search.genome)
    }
    normal.xome <- seqlevels(granges.object)[(regexpr("_", seqlevels(granges.object)) < 0)]
    #xome.value <- str_extract(normal.xome, "[0-9]|1[0-9]|2[0-9]|3[0-9]|4[0-9]|5[0-9]|6[0-9]|7[0-9]|8[0-9]|9[0-9]|X|Y|M")
    positions <- sapply(paste0(normal.xome, "$"), grep, seqnames(seqinfo(search.genome)))
    new2oldmap <- rep(NA, length(seqinfo(search.genome)))
    new2oldmap[positions] <- 1:length(positions)
    seqinfo(granges.object, new2old = new2oldmap) <- seqinfo(search.genome)
  }
  return(granges.object)
}


strSort <- function(x) {
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse = "")
}

#' @importFrom rtracklayer import.bed
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @importFrom GenomicRanges findOverlaps queryHits subjectHits
#' @importFrom GenomeInfoDb species sortSeqlevels
#' @importFrom BiocGenerics sapply
#' @export
snps.from.bed <- function(bedfile = NULL, dbSNP = NULL, search.genome = NULL) {
  snps <- import.bed(bedfile)
  if (Reduce("|", grepl("rs", snps$name)) & (class(dbSNP) != "SNPlocs")) {
    stop(paste0(bedfile, " contains at least one variant with an rsID and no SNPlocs has been indicated\n",
                "Please run availible.SNPs() to check for availble SNPlocs"))
  }
  if (class(search.genome) != "BSgenome") {
    stop(paste0(search.genome, " is not a BSgenome object.\n", "Run availible.genomes() and choose the appropriate BSgenome object"))
  }
  ## spit snps into named and unnamed snps
  snps.noid <- snps[!grepl("rs", snps$name), ]
  ## get ref for unnamed snps
  snps.noid.ref <- getSeq(search.genome, snps.noid)
  snps.noid.ref <- as.character(snps.noid.ref)
  ## get alt for unnamed snps
  snps.noid.alt <- snps.noid$name
  snps.noid.alt <- unlist(lapply(snps.noid.alt, strsplit, split = ":"), recursive = FALSE)
  snps.noid.alt <- sapply(snps.noid.alt, "[", 3)
  ## check if alt was given for unnamed snps
  alt.allele.is.valid <- (toupper(snps.noid.alt) %in% c("A", "T", "G", "C")) &
    (snps.noid.alt != snps.noid.ref)
  if (!Reduce("&", alt.allele.is.valid)) {
    snpnames <- snps.noid$name[!(toupper(snps.noid.alt) %in% c("A", "T", "G",
                                                               "C"))]
    if (length(snpnames) < 50 && length(snpnames) > 0) {
      warning(paste("User variant", snpnames, "alternate allele is not one of \"A\", \"T\", \"G\", or \"C\""))
    } else {
      if (length(snpnames) > 0) {
        warning(paste0(length(snpnames), " user variants contain an alternate allele that is not one of \"A\", \"T\", \"G\", \"C\"\n",
                       " These variants were excluded"))
      }
    }
    equal.to.ref <- snps.noid.alt == snps.noid.ref
    if (sum(equal.to.ref) > 0) {
      warning(paste0(sum(equal.to.ref), " user variants are the same as the reference genome ",
                     search.genome@provider_version, " for ", search.genome@species, "\n These variants were excluded"))
    }
    snps.noid <- snps.noid[alt.allele.is.valid]
    snps.noid.ref <- snps.noid.ref[alt.allele.is.valid]
    snps.noid.alt <- snps.noid.alt[alt.allele.is.valid]
  }
  ## if alt was given calculate ambiguous base
  snps.noid.ambi <- strSort(paste0(snps.noid.alt, snps.noid.ref))
  snps.noid.ambi <- names(IUPAC_CODE_MAP[sapply(as.list(sapply(snps.noid.ambi,
                                                               grep, IUPAC_CODE_MAP)), "[", 1)])
  ## are unnamed snps found in dbsnp ?
  if (species(search.genome) == "Human" && class(dbSNP) == "SNPlocs") {
    snps.noid.chrom <- as.character(seqnames(snps.noid))
    snps.noid.chrom <- unique(snps.noid.chrom)
    snps.noid.chrom <- gsub("chr", "ch", snps.noid.chrom)
    all.dbsnp.chrom <- snplocs(dbSNP, snps.noid.chrom, as.GRanges = T)
    all.dbsnp.chrom <- change.to.search.genome(all.dbsnp.chrom, search.genome)
    present.in.dbsnp <- findOverlaps(snps.noid, all.dbsnp.chrom)
    ## In dbsnp
    dbsnp.for.noid <- all.dbsnp.chrom[subjectHits(present.in.dbsnp), ]
    matches.dbsnp <- snps.noid.ambi[queryHits(present.in.dbsnp)] == dbsnp.for.noid$alleles_as_ambig
    copy.dbsnp <- queryHits(present.in.dbsnp)[matches.dbsnp]
    ## matches dbsnp
    no.dbsnp <- snps.noid[-copy.dbsnp]
    no.dbsnp <- change.to.search.genome(no.dbsnp, search.genome)
  } else {
    copy.dbsnp <- -(1:length(snps.noid))
    matches.dbsnp <- FALSE
    no.dbsnp <- sortSeqlevels(snps.noid)
    no.dbsnp <- change.to.search.genome(no.dbsnp, search.genome)
  }
  mcols(no.dbsnp) <- mcols(no.dbsnp)$name
  colnames(mcols(no.dbsnp)) <- "SNP_id"
  no.dbsnp$alleles_as_ambig <- snps.noid.ambi[-copy.dbsnp]
  no.dbsnp$REF <- snps.noid.ref[-copy.dbsnp]
  no.dbsnp$ALT <- snps.noid.alt[-copy.dbsnp]
  names(no.dbsnp) <- no.dbsnp$SNP_id
  if (Reduce("|", matches.dbsnp)) {
    dbsnp.for.noid <- dbsnp.for.noid[copy.dbsnp, ]
    dbsnp.for.noid$REF <- snps.noid.ref[copy.dbsnp]
    dbsnp.for.noid$ALT <- snps.noid.alt[copy.dbsnp]
    colnames(mcols(dbsnp.for.noid))[1] <- "SNP_id"
    names(dbsnp.for.noid) <- gsub("^", "rs", dbsnp.for.noid$SNP_id)
    warning(paste0("rs", dbsnp.for.noid$SNP_id, " was found as a match for ",
                   snps.noid$name[copy.dbsnp], "; using entry from dbSNP"))
    no.dbsnp <- c(dbsnp.for.noid, no.dbsnp)
  }
  no.dbsnp$REF <- DNAStringSet(no.dbsnp$REF)
  no.dbsnp$ALT <- DNAStringSet(no.dbsnp$ALT)
  no.dbsnp$alleles_as_ambig <- DNAStringSet(no.dbsnp$alleles_as_ambig)
  snps.rsid <- snps[grepl("rs", snps$name), ]
  ## get object for named snps
  if (length(snps.rsid) > 0) {
    snps.rsid.out <- snps.from.rsid(snps.rsid$name, dbSNP = dbSNP, search.genome = search.genome)
    colnames(mcols(snps.rsid.out))[1] <- "SNP_id"
    names(snps.rsid.out) <- gsub("^", "rs", snps.rsid.out$SNP_id)
    snps.out <- c(snps.rsid.out, no.dbsnp)
  } else {
    snps.out <- no.dbsnp
  }
  attributes(snps.out)$genome.package <- attributes(search.genome)$pkgname
  return(snps.out)
}
