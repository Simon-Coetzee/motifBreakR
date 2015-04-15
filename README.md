## motifbreakR
-----

##### Abstract
Functional annotation represents a key step toward the understanding and
interpretation of germline and somatic variation as revealed by genome wide
association studies (GWAS) and The Cancer Genome Atlas (TCGA), respectively.
GWAS have revealed numerous genetic risk variants residing in non-coding DNA
associated with complex diseases. For sequences that lie within enhancers or
promoters of transcription, it is straightforward to assess the effects of
variants on likely transcription factor binding sites. We introduce
motifbreakR, which allows the biologist to judge whether the sequence
surrounding a polymorphism or mutation is a good match, and how much
information is gained or lost in one allele of the polymorphism relative to
another or mutation vs. wildtype. MotifbreakR is flexible, giving a choice of
algorithms for interrogation of genomes with motifs from public sources that
users can choose from; these are 1) a weighted-sum, 2) log-probabilities, and
3) relative entropy. MotifbreakR can predict effects for novel or previously
described variants in public databases, making it suitable for tasks beyond
the scope of its original design. Lastly, it can be used to interrogate any
genome curated within bioconductor.

### Install

```{r}
# Install prerequisite packages from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(c("BiocParallel", "motifStack", "BSgenome", "BiocGenerics",
           "Biostrings", "GenomeInfoDb", "GenomicRanges", "Gviz", "S4Vectors",
           "rtracklayer", "IRanges", "BSgenome.Hsapiens.UCSC.hg19",
           "SNPlocs.Hsapiens.dbSNP.20120608"))
# Install motifbreakR from github
install.packages("devtools")
devtools::install_github("Simon-Coetzee/MotifBreakR")
```

Thus far the vignette is incomplete but see the help for `motifbreakR` and 
`plotMB` to get started.
