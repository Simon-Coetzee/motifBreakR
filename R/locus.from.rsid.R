#' Import SNPs from rsid for use in motifbreakR
#'
#' @param rsid Character; a character vector of rsid values from dbSNP
#' @param dbSNP an object of class SNPlocs to lookup rsids; see \code{availible.SNPs} in
#'   \code{\link[BSgenome]{injectSNPs}} to check for availible SNPlocs
#' @param search.genome an object of class BSgenome for the species you are interrogating;
#'  see \code{\link[BSgenome]{available.genomes}} for a list of species
#' @seealso See \code{\link{motifbreakR}} for analysis; See \code{\link{snps.from.file}}
#'   for an alternate method for generating a list of variants.
#' @details \code{snps.from.rsid} take an rsid, or character vector of rsids and
#'  generates the required object to input into \code{motifbreakR}
#' @return a GRanges object containing:
#'  \item{SNP_id}{The rsid of the snp with the "rs" portion stripped}
#'  \item{alleles_as_ambig}{THE IUPAC ambiguity code between the reference and
#'  alternate allele for this SNP}
#'  \item{REF}{The reference allele for the SNP}
#'  \item{ALT}{The alternate allele for the SNP}
#'  @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  library(SNPlocs.Hsapiens.dbSNP.20120608)
#'  snps.file <- system.file("extdata", "pca.enhancer.snps", package = "motifbreakR")
#'  snps <- as.character(read.table(snps.file)[,1])
#'  snps.mb <- snps.from.rsid(snps,
#'                            dbSNP = SNPlocs.Hsapiens.dbSNP.20120608,
#'                            search.genome = BSgenome.Hsapiens.UCSC.hg19)
#'
#' @importFrom BSgenome snpid2grange snplocs
#' @importFrom Biostrings DNAStringSet
#' @export
snps.from.rsid <- function(rsid = NULL, dbSNP = NULL,
                           search.genome = NULL) {
  if (is.null(rsid)) {
    stop("no RefSNP IDs have been found, please include RefSNP ID numbers")
  }
  if (!inherits(dbSNP, "SNPlocs")) {
    stop(paste0("dbSNP argument was not provided with a valid SNPlocs object.\n",
                "Please run availible.SNPs() to check for availible SNPlocs"))
  }
  if (!inherits(search.genome, "BSgenome")) {
    stop(paste0("search.genome argument was not provided with a valid BSgenome object.\n",
                "Run availible.genomes() and choose the appropriate BSgenome object"))
  }
  if (Reduce("&", !grepl("rs", rsid))) {
    bad.names <- rsid[!grepl("rs", rsid)]
    stop(paste(paste(bad.names, collapse = " "), "are not rsids, perhaps you want to import your snps from a bed or vcf file with snps.from.file()?"))
  }
  rsid <- unique(rsid)
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

#' @import GenomeInfoDb
#' @importFrom stringr str_extract
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
    positions <- unlist(sapply(paste0(normal.xome, "$"), grep, seqnames(seqinfo(search.genome))))
    new2oldmap <- rep(NA, length(seqinfo(search.genome)))
    new2oldmap[positions] <- 1:length(positions)
    seqinfo(granges.object, new2old = new2oldmap) <- seqinfo(search.genome)
  }
  return(granges.object)
}


strSort <- function(x) {
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse = "")
}
#' Import SNPs from a BED file or VCF file for use in motifbreakR
#'
#' @param file Character; a character containing the path to a bed file or a vcf file
#'   see Details for a description of the required format
#' @param dbSNP OPTIONAL; an object of class SNPlocs to lookup rsids; see \code{availible.SNPs} in
#'   \code{\link[BSgenome]{injectSNPs}} to check for availible SNPlocs
#' @param search.genome an object of class BSgenome for the species you are interrogating;
#'  see \code{\link[BSgenome]{available.genomes}} for a list of species
#' @param format Character; one of \code{bed} or \code{vcf}
#' @seealso See \code{\link{motifbreakR}} for analysis; See \code{\link{snps.from.rsid}}
#'   for an alternate method for generating a list of variants.
#' @details \code{snps.from.file} takes a character vector describing the file path
#'  to a bed file that contains the necissary information to generate the input for
#'  \code{motifbreakR} see \url{http://www.genome.ucsc.edu/FAQ/FAQformat.html#format1}
#'  for a complete description of the BED format.  Our convention deviates in that there
#'  is a required format for the name field.  \code{name} is defined as chromosome:start:REF:ALT
#'  or the rsid from dbSNP (if you've included the optional SNPlocs argument).
#'  For example if you were to include rs123 in it's alternate
#'  format it would be entered as chr7:24966446:C:A
#' @return a GRanges object containing:
#'  \item{SNP_id}{The rsid of the snp with the "rs" portion stripped}
#'  \item{alleles_as_ambig}{THE IUPAC ambiguity code between the reference and
#'  alternate allele for this SNP}
#'  \item{REF}{The reference allele for the SNP}
#'  \item{ALT}{The alternate allele for the SNP}
#'  @examples
#'  library(BSgenome.Drerio.UCSC.danRer7)
#'  library(SNPlocs.Hsapiens.dbSNP.20120608)
#'  snps.bed.file <- system.file("extdata", "danRer.bed", package = "motifbreakR")
#'  # see the contents
#'  read.table(snps.bed.file, header = FALSE)
#'  #import the BED file
#'  snps.mb <- snps.from.file(snps.bed.file,
#'                            search.genome = BSgenome.Drerio.UCSC.danRer7,
#'                            format = "bed")
#'
#' @importFrom rtracklayer import
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @importFrom VariantAnnotation readVcf
#' @export
snps.from.file <- function(file = NULL, dbSNP = NULL, search.genome = NULL, format = "bed") {
  if(format == "vcf"){
    if (!inherits(search.genome, "BSgenome")) {
      stop(paste0(search.genome, " is not a BSgenome object.\n", "Run availible.genomes() and choose the appropriate BSgenome object"))
    }
    genome.name <- genome(search.genome)[[1]]
    snps <- rowRanges(readVcf(file, genome.name))
    ALTS <- unlist(snps$ALT)
    numsplits <- elementLengths(snps$ALT)
    snps <- rep(snps, times = numsplits)
    snps$ALT <- ALTS
    snps$ALT <- as.character(snps$ALT)
    snps$REF <- as.character(snps$REF)
    snps.ambi <- strSort(paste0(snps$ALT, snps$REF))
    IUPAC_code_revmap <- names(IUPAC_CODE_MAP)
    names(IUPAC_code_revmap) <- IUPAC_CODE_MAP
    snps.ambi <- IUPAC_code_revmap[snps.ambi]
    names(snps.ambi) <- NULL
    snps$alleles_as_ambig <- snps.ambi
    snps$SNP_id <- names(snps)
    mcols(snps) <- mcols(snps)[, c("SNP_id", "alleles_as_ambig", "REF", "ALT")]
    snps$REF <- DNAStringSet(snps$REF)
    snps$ALT <- DNAStringSet(snps$ALT)
    snps$alleles_as_ambig[is.na(snps$alleles_as_ambig)] <- ""
    snps$alleles_as_ambig <- DNAStringSet(snps$alleles_as_ambig)
    snps <- snps[seqnames(snps) != "MT"]
    seqlevels(snps) <- paste0("chr", seqlevels(snps))
    snps <- change.to.search.genome(snps, search.genome)
    can.ref <- getSeq(search.genome, snps)
    names(can.ref) <- NULL
    if(isTRUE(all.equal(can.ref, snps$REF))) {
      rm(can.ref)
    } else {
      warning(paste0("User selected reference allele differs from the sequence in ",
                     attributes(search.genome)$pkgname, " continuing with genome specified",
                     " reference allels\n", "there are ", sum(snps$REF != can.ref),
                     " differences"))
    }
    attributes(snps)$genome.package <- attributes(search.genome)$pkgname
    return(snps)
  } else {
    if(format == "bed") {
      snps <- import(file, format = "bed")
      if (Reduce("|", grepl("rs", snps$name)) & (!inherits(dbSNP, "SNPlocs"))) {
        stop(paste0(file, " contains at least one variant with an rsID and no SNPlocs has been indicated\n",
                    "Please run availible.SNPs() to check for availble SNPlocs"))
      }
      if (!inherits(search.genome, "BSgenome")) {
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
      snps.noid.ref.user <- sapply(snps.noid.alt, "[", 3)
      if(isTRUE(all.equal(snps.noid.ref, snps.noid.ref.user))) {
        rm(snps.noid.ref.user)
      } else {
        warning(paste0("User selected reference allele differs from the sequence in ",
                       attributes(search.genome)$pkgname, " continuing with genome specified",
                       " reference allels\n", " there are ", sum(snps.noid.ref != snps.noid.ref.user),
                       " differences"))
      }
      snps.noid.alt <- sapply(snps.noid.alt, "[", 4)
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
                         search.genome@provider_version, " for ", search.genome@common_name, "\n These variants were excluded"))
        }
        snps.noid <- snps.noid[alt.allele.is.valid]
        snps.noid.ref <- snps.noid.ref[alt.allele.is.valid]
        snps.noid.alt <- snps.noid.alt[alt.allele.is.valid]
      }
      ## if alt was given calculate ambiguous base
      snps.noid.ambi <- strSort(paste0(snps.noid.alt, snps.noid.ref))
      IUPAC_code_revmap <- names(IUPAC_CODE_MAP)
      names(IUPAC_code_revmap) <- IUPAC_CODE_MAP
      snps.noid.ambi <- IUPAC_code_revmap[snps.noid.ambi]
      names(snps.noid.ambi) <- NULL
      #snps.noid.ambi <- names(IUPAC_CODE_MAP[sapply(as.list(sapply(snps.noid.ambi,
      #                                                             grep, IUPAC_CODE_MAP)), "[", 1)])
      ## are unnamed snps found in dbsnp ?
      if (commonName(search.genome) == "Human" && class(dbSNP) == "SNPlocs") {
        snps.noid.chrom <- as.character(seqnames(snps.noid))
        snps.noid.chrom <- unique(snps.noid.chrom)
        snps.noid.chrom <- gsub("chr", "ch", snps.noid.chrom)
        all.dbsnp.chrom <- snplocs(dbSNP, snps.noid.chrom, as.GRanges = TRUE)
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
    } else {
      stop("format must be one of 'vcf' or 'bed'; currently set as ", format)
    }
  }
}

