## motifbreakR
-----

##### Documentation
See the [motifbreakR vignette](https://bioconductor.org/packages/motifbreakR) on Bioconductor, or `vignette("motifbreakR-vignette", package = "motifbreakR")`, for an introduction to `motifbreakR`

See `help("motifbreakR")` for detailed help with running `motifbreakR`.

See `help("plotMB")` for detailed help with visualization.

Please cite:

- Coetzee SG, Hazelett DJ (2024). motifbreakR v2: expanded variant analysis including indels and integrated evidence from transcription factor binding databases. *Bioinformatics Advances*, 4(1), vbae162. [doi:10.1093/bioadv/vbae162](https://doi.org/10.1093/bioadv/vbae162)
- Coetzee SG, Coetzee GA, Hazelett DJ (2015). motifbreakR: an R/Bioconductor package for predicting variant effects at transcription factor binding sites. *Bioinformatics*, 31(23), 3847-3849. [doi:10.1093/bioinformatics/btv470](https://doi.org/10.1093/bioinformatics/btv470)

##### Abstract
Functional annotation represents a key step toward the understanding and
interpretation of germline and somatic variation as revealed by genome wide
association studies (GWAS) and The Cancer Genome Atlas (TCGA), respectively.
GWAS have revealed numerous genetic risk variants residing in non-coding DNA
associated with complex diseases. For sequences that lie within enhancers or
promoters of transcription, it is straightforward to assess the effects of
variants on likely transcription factor binding sites. We introduce
_motifbreakR_, which allows the biologist to judge whether the sequence
surrounding a polymorphism or mutation (a single nucleotide variant or a small
insertion or deletion) is a good match, and how much
information is gained or lost in one allele of the polymorphism relative to
another or mutation vs. wildtype. _MotifbreakR_ is flexible, giving a choice of
algorithms for interrogation of genomes with motifs from public sources that
users can choose from; these are 1) a weighted-sum, 2) log-probabilities, and
3) relative entropy. _MotifbreakR_ can predict effects for novel or previously
described variants in public databases, making it suitable for tasks beyond
the scope of its original design. Lastly, it can be used to interrogate any
genome curated within Bioconductor.

### Install

_motifbreakR_ is available from [Bioconductor](https://bioconductor.org/packages/motifbreakR):
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("motifbreakR")
```

The examples and vignette also use genome and SNP annotation packages, for example:
```r
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19",
                       "SNPlocs.Hsapiens.dbSNP155.GRCh37"))
```

#### Development version
The development version can be installed from GitHub:
```r
BiocManager::install("Simon-Coetzee/motifBreakR")
```
