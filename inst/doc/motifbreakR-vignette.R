## ----style, echo = FALSE, results='hide', message=FALSE, warning=FALSE----------------------------
BiocStyle::markdown()
options(width=100)
knitr::opts_chunk$set(cache=TRUE, autodep=TRUE)

## ----start, cache=TRUE, message=FALSE-------------------------------------------------------------
library(motifbreakR)
pca.snps.file <- system.file("extdata", "pca.enhancer.snps", package = "motifbreakR")
pca.snps <- as.character(read.table(pca.snps.file)[,1])

## ----outline, eval=FALSE--------------------------------------------------------------------------
#  variants <- snps.from.rsid(rsid = pca.snps,
#                             dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
#                             search.genome = BSgenome.Hsapiens.UCSC.hg19)
#  motifbreakr.results <- motifbreakR(snpList = variants, pwmList = MotifDb, threshold = 0.9)
#  plotMB(results = motifbreakr.results, rsid = "rs7837328", effect = "strong")

## ----whichsnps, message=FALSE, cache=TRUE---------------------------------------------------------
library(BSgenome)
available.SNPs()

## ----showbed, message=FALSE, cache=TRUE-----------------------------------------------------------
snps.file <- system.file("extdata", "snps.bed", package = "motifbreakR")
read.delim(snps.file, header = FALSE)

## ----fromrsid, message=FALSE, cache=TRUE----------------------------------------------------------
library(SNPlocs.Hsapiens.dbSNP142.GRCh37) # dbSNP142 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)     # hg19 genome
head(pca.snps)
snps.mb <- snps.from.rsid(rsid = pca.snps,
                          dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)
snps.mb

## ----whichgenome, message=FALSE, cache=TRUE-------------------------------------------------------
library(BSgenome)
genomes <- available.genomes()
length(genomes)
genomes

## ----showbed2, message=TRUE, cache=TRUE-----------------------------------------------------------
snps.bed.file <- system.file("extdata", "snps.bed", package = "motifbreakR")
# see the contents
read.table(snps.bed.file, header = FALSE)

## ----getbedrs, message=TRUE, cache=TRUE-----------------------------------------------------------
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
#import the BED file
snps.mb.frombed <- snps.from.file(file = snps.bed.file,
                                  dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                  format = "bed")
snps.mb.frombed

## ----getbedcust, message=FALSE, cache=TRUE--------------------------------------------------------
library(BSgenome.Drerio.UCSC.danRer7)
snps.bedfile.nors <- system.file("extdata", "danRer.bed", package = "motifbreakR")
read.table(snps.bedfile.nors, header = FALSE)
snps.mb.frombed <- snps.from.file(file = snps.bedfile.nors,
                                  search.genome = BSgenome.Drerio.UCSC.danRer7,
                                  format = "bed")
snps.mb.frombed

## ----motifdb, message=FALSE, cache=TRUE-----------------------------------------------------------
library(MotifDb)
MotifDb

## ----motifdbtableshow, message=FALSE, cache=TRUE, eval=FALSE--------------------------------------
#  ### Here we can see which organisms are availible under which sources
#  ### in MotifDb
#  table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource)

## ----motifdbtableres, message=FALSE, cache=TRUE, echo=FALSE---------------------------------------
knitr::kable(table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource), format = "html", table.attr="class=\"table table-striped table-hover\"")

## ----motifbreakrmot, message=FALSE, cache=TRUE----------------------------------------------------
data(motifbreakR_motif)
motifbreakR_motif

## ----hocomocomot, message=FALSE, cache=TRUE-------------------------------------------------------
data(hocomoco)
hocomoco

## ----runmotifbreakr, message=TRUE, cache=TRUE-----------------------------------------------------
results <- motifbreakR(snpList = snps.mb[1:5], filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())

## ----getonesnp, cache=TRUE------------------------------------------------------------------------
rs1006140 <- results[names(results) %in% "rs1006140"]
rs1006140

## ----calcp.show, cache=TRUE, eval=FALSE-----------------------------------------------------------
#  rs1006140 <- calculatePvalue(rs1006140)
#  rs1006140

## ----calcp, cache=TRUE, echo=FALSE----------------------------------------------------------------
data(example.pvalue.rda)
rs1006140

## ----aboutparallel, message=TRUE, cache=TRUE------------------------------------------------------
BiocParallel::registered()
BiocParallel::bpparam()

## ----plotting, cache=TRUE, fig.retina=2, fig.align='center', fig.height=8, fig.width=6------------
plotMB(results = results, rsid = "rs1006140", effect = "strong")

## ---- cache=TRUE----------------------------------------------------------------------------------
sessionInfo()

