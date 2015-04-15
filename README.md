## motifbreakR
-----

### Install
```{r}
install.packages("devtools")
source("http://bioconductor.org/biocLite.R")
biocLite(c("BiocParallel", "motifStack", "BSgenome", "BiocGenerics",
           "Biostrings", "GenomeInfoDb", "GenomicRanges", "Gviz", "S4Vectors",
           "rtracklayer", "IRanges", "BSgenome.Hsapiens.UCSC.hg19",
           "SNPlocs.Hsapiens.dbSNP.20120608"))
devtools::install_github("Simon-Coetzee/MotifBreakR")
```

Thus far the vignette is incomplete but see the help for `motifbreakR` and 
`plotMB` to get started.
