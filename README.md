## motifbreakR
-----

##### Documentation
See [motifbreakR vignette](http://simon-coetzee.github.io/motifBreakR/) for an introduction to `motifbreakR`

See `help("motifbreakR")` for detailed help with running `motifbreakR`.

See `help("plotMB")` for detailed help with visualization.

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
genome curated within Bioconductor.

### Install
#### Prepairing to install

You will need ghostscript: the full path to the executable can be set by the environment variable R_GSCMD. If this is unset a GhostScript executable will be searched by name on your path. For example on a Unix, Linux, or Mac environment `gs` will be used for searching, and on a Windows environment the setting of the environment variable GSC is used, otherwise commands `gswi64c.exe` and then `gswin32c.exe` are tried.

For example on Windows, assume that the gswin32c.exe is installed at `C:\Program Files\gs\gs9.06\bin`, open R then try:
```{r}
Sys.setenv(R_GSCMD="\"C:\\Program Files\\gs\\gs9.06\\bin\\gswin32c.exe\"")
```

In Linux try your package manager to search for ghostscript.
In Mac [Homebrew](http://brew.sh/) serves as a great package manager.
Windows can find it [here](http://ghostscript.com/download/gsdnld.html)

Additionally, one of the packages that we depend on depends upon [MotIV](http://www.bioconductor.org/packages/release/bioc/html/MotIV.html) in Bioconductor which in turn depends upon the GNU Scientific Library.
Please see the MotIV Vignette, appendix: [GSL Installation](http://www.bioconductor.org/packages/release/bioc/vignettes/MotIV/inst/doc/MotIV.pdf#section.11) If you are having issues with installation.

If you'd like to make the vignette yourself (and have it appear identically to the one listed above), you also need [pandoc](http://pandoc.org/installing.html) which converts R markdown into any of numerous formats including `.html`, `.pdf`, and microsoft word style `.doc` files.

#### Getting prerequisite packages from Bioconductor
```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(c("BiocParallel", "motifStack", "BSgenome", "BiocGenerics",
           "Biostrings", "GenomeInfoDb", "GenomicRanges", "Gviz", "S4Vectors",
           "rtracklayer", "IRanges", "MotifDb", "BSgenome.Hsapiens.UCSC.hg19",
           "SNPlocs.Hsapiens.dbSNP.20120608", "SNPlocs.Hsapiens.dbSNP155.GRCh37",
           "VariantAnnotation", "matrixStats", "BiocStyle"))
install.packages(c("TFMPvalue", "knitr", "rmarkdown"))
```

#### Install motifbreakR from github
```{r}
install.packages("devtools")
devtools::install_github("Simon-Coetzee/motifBreakR")
```
