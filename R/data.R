#' MotifDb object containing motif information from Homo Sapiens Comprehensive
#' Model Collection (HOCOMOCO) of transcription factor (TF) binding models
#'
#' From the abstract:
#' "We present the Homo sapiens comprehensive model collection (HOCOMOCO,
#' \url{http://autosome.ru/HOCOMOCO/}, \url{http://cbrc.kaust.edu.sa/hocomoco/})
#' containing carefully hand-curated TFBS models constructed by integration of
#' binding sequences obtained by both low- and high-throughput methods. To
#' construct position weight matrices to represent these TFBS models, we used
#' ChIPMunk software in four computational modes, including newly developed
#' periodic positional prior mode associated with DNA helix pitch. We selected
#' only one TFBS model per TF, unless there was a clear experimental evidence
#' for two rather distinct TFBS models. We assigned a quality rating to each
#' model. HOCOMOCO contains 426 systematically curated TFBS models for 401 human
#' TFs, where 172 models are based on more than one data source."
#'
#' Load with \code{data(hocomoco)}
#'
#'@return \code{\link[MotifDb]{MotifList-class}} object
#'
#'@format \code{\link[MotifDb]{MotifDb}} object of length 426; to access metadata
#' use mcols(hocomoco)
#'\describe{
#'  \item{providerName}{Name provided by HOCOMOCO}
#'  \item{providerId}{ID provided by HOCOMOCO including experiment type}
#'  \item{dataSource}{\code{"HOCOMOCO"}}
#'  \item{geneSymbol}{Gene symbol for the transcription factor}
#'  \item{geneId}{Entrez gene id for the transcription factor}
#'  \item{geneIdType}{\code{"ENTREZ"}}
#'  \item{proteinId}{UNIPROT id for the transcription factor}
#'  \item{proteinIdType}{\code{"UNIPROT"}}
#'  \item{organism}{\code{"Hsapiens"}}
#'  \item{sequenceCount}{Number of sequences evaluated for producing the PWM}
#'  \item{bindingSequence}{Consensus sequence for the motif}
#'  \item{bindingDomain}{\code{NA} incomplete}
#'  \item{tfFamily}{\code{NA} incomplete}
#'  \item{experimentType}{from \url{http://autosome.ru/HOCOMOCO/Details.php#200}
#'    quoted here:
#'
#' "TFBS model identification modes
#'
#' To construct TFBS models ChIPMunk was run four times: two times (f1) and (f2)
#' with uniform model positional prior and two times (si) and (do) with
#' informative model positional prior.
#'
#' The min-to-max (f1) model length estimation mode was used with the min length
#' of 7 bp and increasing it by 1 bp until the default max length of 25 bp was
#' reached following the optimal length selection procedure as in Kulakovskiy
#' and Makeev, Biophysics, 2009. For max-to-min (f2) model length estimation
#' mode we started from 25 bp and searched for the best alignment decreasing the
#' length by 1 bp until the minimal length of 7 bp. We also used the single (si)
#' and double box (do) model positional priors in order to simulate DNA helix
#' turn. For a single box, the positional weights are to be distributed as
#' cos2(\eqn{pi} n / T), where T=10.5 is the DNA helix pitch, n is the coordinate
#' within the alignment, and the center of the alignment of the length L is at
#' n=0. During the internal cycle of PWM optimization the PWM column scores are
#' multiplied by prior values so the columns closer to the center of the
#' alignment (n=0) receive no score penalty while the columns around (n =
#' 5,6,-5,-6) contribute much less to the score of the PWM under optimization.
#' The single box model prior was used along with the min-to-max length
#' estimation mode (si). We also used the double box model prior with a shape
#' prior equal to sin2(\eqn{pi}n / T), which was used to search for possibly longer
#' double box models in the max-to-min length estimation mode (do).
#'
#' Model quality assignment
#'
#' The resulting models were rated (from A to F) according to their quality.
#' Model quality rates from A-to-D were assigned to proteins known to be TFs,
#' including those listed in Schaefer et al., Nucleic Acids Research, 2011 with
#' addition of a number of proteins having relevant models and sufficient
#' evidence to be TFs. The ratings were assigned by human curation according to
#' the following criteria:
#'
#' Relevant distribution of position-specific information content over alignment
#' columns, which means a model LOGO representation displaying well formed core
#' positions with a high information content surrounded by flanking letters with
#' lower information content; the information content at flanking positions
#' decreasing with the distance from the model core.
#'
#' "Stability", which means that in more than one of the ChIPMunk modes we
#' obtained models with a similar length, consensus, and comparable number of
#' aligned binding sites, along with a similar shape of model LOGO
#' representation. "Similarity" of the model to the binding sequence consensus
#' for this TF given in the UniProt or other databases, which means similarity
#' of the shape of the model LOGO and TFBS lengths to those of other TFs from
#' the same TF family. "A total number of binding sites" was also considered as
#' a quality measure, as a large set of binding regions (mostly but not limited
#' to ChIP-Seq and parallel SELEX) implies that there are many observations of
#' each letter in any position of the alignment, particularly many observations
#' of non-consensus letters in core positions. In positions with low information
#' content, where there is no strong consensus, all variants have many
#' observations, and thus the observed letter frequencies are less dependent on
#' statistical fluctuations.
#'
#' Quality A was assigned to high confidence models complying with all four
#' criteria listed in the section above. Quality B was assigned to models built
#' from large sequence sets that failed no more than one out of the three
#' remaining criteria. Quality C was assigned to models built from small
#' sequence sets but (with a number of specifically marked exceptions) complying
#' with the three remaining criteria. Quality D models missed part of the known
#' consensus sequence or had no clearly significant core positions in the TFBS
#' model. Quality E (error) was assigned to models for proteins not convincingly
#' shown to be TFs or to models exhibiting an irrelevant LOGO shape or a wrong
#' consensus sequence. Quality F (failure) was assigned to TFs for which there
#' was no reliable model identified."
#'}
#'  \item{pubmedID}{\code{"23175603"} see \code{Source} for more details}
#'}
#'
#' @source Kulakovskiy,I.V., Medvedeva,Y.A., Schaefer,U., Kasianov,A.S.,
#'   Vorontsov,I.E., Bajic,V.B. and Makeev,V.J. (2013) HOCOMOCO: a comprehensive
#'   collection of human transcription factor binding sites models. Nucleic
#'   Acids Research, \bold{41}, D195--D202.
#'
#' @seealso \url{http://autosome.ru/HOCOMOCO/} \url{http://cbrc.kaust.edu.sa/hocomoco/}
#' @examples
#' data(hocomoco)
#' hocomoco
"hocomoco"

#' MotifDb object containing motif information from the known and discovered
#' motifs for the ENCODE TF ChIP-seq datasets.
#'
#' From the abstract: "Recent advances in technology have led to a dramatic
#' increase in the number of available transcription factor ChIP-seq and
#' ChIP-chip data sets. Understanding the motif content of these data sets is an
#' important step in understanding the underlying mechanisms of regulation. Here
#' we provide a systematic motif analysis for 427 human ChIP-seq data sets using
#' motifs curated from the literature and also discovered de novo using five
#' established motif discovery tools. We use a systematic pipeline for
#' calculating motif enrichment in each data set, providing a principled way for
#' choosing between motif variants found in the literature and for flagging
#' potentially problematic data sets. Our analysis confirms the known
#' specificity of 41 of the 56 analyzed factor groups and reveals motifs of
#' potential cofactors. We also use cell type-specific binding to find factors
#' active in specific conditions. The resource we provide is accessible both for
#' browsing a small number of factors and for performing large-scale systematic
#' analyses. We provide motif matrices, instances and enrichments in each of the
#' ENCODE data sets. The motifs discovered here have been used in parallel
#' studies to validate the specificity of antibodies, understand cooperativity
#' between data sets and measure the variation of motif binding across
#' individuals and species."
#'
#' Load with \code{data(encodemotif)}
#'
#'@return \code{\link[MotifDb]{MotifList-class}} object
#'
#'@format \code{\link[MotifDb]{MotifDb}} object of length 2064; to access metadata
#' use mcols(encodemotif)
#'\describe{
#'  \item{providerName}{Name provided by ENCODE}
#'  \item{providerId}{Same as providerName}
#'  \item{dataSource}{\code{"ENCODE-motif"}}
#'  \item{geneSymbol}{Gene symbol for the transcription factor}
#'  \item{geneId}{Entrez gene id for the transcription factor}
#'  \item{geneIdType}{\code{"ENTREZ"}}
#'  \item{proteinId}{UNIPROT id for the transcription factor}
#'  \item{proteinIdType}{\code{"UNIPROT"}}
#'  \item{organism}{\code{"Hsapiens"}}
#'  \item{sequenceCount}{\code{NA} not available}
#'  \item{bindingSequence}{Consensus sequence for the motif}
#'  \item{bindingDomain}{\code{NA} incomplete}
#'  \item{tfFamily}{\code{NA} incomplete}
#'  \item{experimentType}{occurs in two forms:
#'
#'  For motifs that were discovered in this study, the format is \code{cellType_source-LabMetadata:MotifFinder#Location} for example \code{H1-hESC_encode-Myers_seq_hsa_v041610.2_r1:MEME#2#Intergenic}.
#'
#'  For motifs that were "known" the format tends to be \code{TF_source_sourceId} for example \code{AP1_jaspar_MA0099.2}.
#'}
#'  \item{pubmedID}{\code{"24335146"} see \code{Source} for more details}
#'}
#'
#' @source Pouya Kheradpour and Manolis Kellis (2013 December 13) Systematic
#'   discovery and characterization of regulatory motifs in ENCODE TF binding
#'   experiments. Nucleic Acids Research, doi:10.1093/nar/gkt1249
#'
#' @seealso \url{http://compbio.mit.edu/encode-motifs/}
#' @examples
#' data(encodemotif)
#' encodemotif
"encodemotif"

#' MotifDb object containing motif information from around the genomic regions
#' bound by 119 human transcription factors in Factorbook.
#'
#' From the abstract: "Chromatin immunoprecipitation coupled with
#' high-throughput sequencing (ChIP-seq) has become the dominant technique for
#' mapping transcription factor (TF) binding regions genome-wide. We performed
#' an integrative analysis centered around 457 ChIP-seq data sets on 119 human
#' TFs generated by the ENCODE Consortium. We identified highly enriched
#' sequence motifs in most data sets, revealing new motifs and validating known
#' ones. The motif sites (TF binding sites) are highly conserved evolutionarily
#' and show distinct footprints upon DNase I digestion. We frequently detected
#' secondary motifs in addition to the canonical motifs of the TFs, indicating
#' tethered binding and cobinding between multiple TFs. We observed significant
#' position and orientation preferences between many cobinding TFs. Genes
#' specifically expressed in a cell line are often associated with a greater
#' occurrence of nearby TF binding in that cell line. We observed
#' cell-line-specific secondary motifs that mediate the binding of the histone
#' deacetylase HDAC2 and the enhancer-binding protein EP300. TF binding sites
#' are located in GC-rich, nucleosome-depleted, and DNase I sensitive regions,
#' flanked by well-positioned nucleosomes, and many of these features show cell
#' type specificity. The GC-richness may be beneficial for regulating TF binding
#' because, when unoccupied by a TF, these regions are occupied by nucleosomes
#' in vivo. We present the results of our analysis in a TF-centric web
#' repository Factorbook (\url{http://factorbook.org}) and will continually update
#' this repository as more ENCODE data are generated."
#'
#' Load with \code{data(factorbook)}
#'
#'@return \code{\link[MotifDb]{MotifList-class}} object
#'
#'
#'@format \code{\link[MotifDb]{MotifDb}} object of length 79; to access metadata
#' use mcols(factorbook)
#'\describe{
#'  \item{providerName}{Name listed in meme output of 'Supp TableS2.pdf' for the citation indicated below}
#'  \item{providerId}{Same as providerName}
#'  \item{dataSource}{\code{"FactorBook"}}
#'  \item{geneSymbol}{\code{NA} these motifs don't have a direct 1 to 1 relationship with a transcription factor}
#'  \item{geneId}{\code{NA}}
#'  \item{geneIdType}{\code{NA}}
#'  \item{proteinId}{\code{NA}}
#'  \item{proteinIdType}{\code{NA}}
#'  \item{organism}{\code{"Hsapiens"}}
#'  \item{sequenceCount}{\code{NA}}
#'  \item{bindingSequence}{Consensus sequence for the motif}
#'  \item{bindingDomain}{\code{NA}}
#'  \item{tfFamily}{\code{NA}}
#'  \item{experimentType}{\code{NA}}
#'  \item{pubmedID}{\code{"22955990"} see \code{Source} for more details}
#'}
#'
#' @source J Wang, J Zhuang, S Iyer, XY Lin, et al. (2012) Sequence features and
#'   chromatin structure around the genomic regions bound by 119 human transcription
#'   factors. Genome Research, \bold{22 (9)}, 1798-1812, doi:10.1101/gr.139105.112
#'
#' @seealso \url{http://factorbook.org}
#' @examples
#' data(factorbook)
#' factorbook
"factorbook"

#' MotifDb object containing motif information from motif databases included in
#' HOMER.
#'
#' From the website: "Homer includes several motif databases that are used to help annotate
#' results and conduct searches for known motifs.  HOMER contains a custom motif
#' database based on independent analysis of mostly ChIP-Seq data sets which is
#' heavily utilized in the software." See \url{http://homer.salk.edu/homer/motif/motifDatabase.html}
#' for more information on how these files were generated, and Homer's sources.
#'
#' Load with \code{data(homer)}
#'
#'@return \code{\link[MotifDb]{MotifList-class}} object
#'
#'
#'@format \code{\link[MotifDb]{MotifDb}} object of length 247; to access metadata
#' use mcols(homer)
#'\describe{
#'  \item{providerName}{Name provided HOMER}
#'  \item{providerId}{Factor Name provided by HOMER}
#'  \item{dataSource}{\code{"HOMER"}}
#'  \item{geneSymbol}{Symbol provided by HOMER}
#'  \item{geneId}{Entrez gene id for the transcription factor}
#'  \item{geneIdType}{\code{"ENTREZ"}}
#'  \item{proteinId}{UNIPROT id for the transcription factor}
#'  \item{proteinIdType}{\code{"UNIPROT"}}
#'  \item{organism}{\code{"Hsapiens"}}
#'  \item{sequenceCount}{\code{NA}}
#'  \item{bindingSequence}{Consensus sequence for the motif}
#'  \item{bindingDomain}{DBD provided by HOMER}
#'  \item{tfFamily}{\code{NA}}
#'  \item{experimentType}{The Celltype, IP, Assay, and GEO id if applicable for the motif}
#'  \item{pubmedID}{\code{"20513432"} see \code{Source} for more details}
#'}
#'
#' @source Heinz S, Benner C, Spann N, Bertolino E et al. (2010 May 28) Simple Combinations
#'   of Lineage-Determining Transcription Factors Prime cis-Regulatory
#'   Elements Required for Macrophage and B Cell Identities. Mol Cell, \bold{38(4):576-589}.
#'   PMID: \href{http://www.ncbi.nlm.nih.gov/sites/entrez?Db=Pubmed&term=20513432[UID]}{20513432}
#'
#' @seealso \url{http://homer.salk.edu/homer/index.html} \url{http://homer.salk.edu/homer/motif/motifDatabase.html}
#'   \url{http://homer.salk.edu/homer/motif/HomerMotifDB/homerResults.html}
#' @examples
#' data(homer)
#' homer
"homer"

#' MotifDb object containing motif information from the motif databases of HOCOMOCO, Homer,
#'  FactorBook and ENCODE
#'
#' This object contains all the \code{\link[MotifDb]{MotifList-class}} objects that were generated
#' for this package.  See the individual help sections for \code{\link{hocomoco}}, \code{\link{homer}},
#' \code{\link{factorbook}}, and \code{\link{encodemotif}}, for how the data is formatted.
#'
#' Load with \code{data(motifbreakR_motif)}
#'
#'@return \code{\link[MotifDb]{MotifList-class}} object
#'
#'@format \code{\link[MotifDb]{MotifDb}} object of length 2816; to access metadata
#' use mcols(motifbreakR_motif)
#'
#' @seealso \code{\link{hocomoco}}, \code{\link{homer}},
#' \code{\link{factorbook}}, and \code{\link{encodemotif}}
#'
#'
#' @source Kulakovskiy,I.V., Medvedeva,Y.A., Schaefer,U., Kasianov,A.S.,
#'   Vorontsov,I.E., Bajic,V.B. and Makeev,V.J. (2013) HOCOMOCO: a comprehensive
#'   collection of human transcription factor binding sites models. Nucleic
#'   Acids Research, \bold{41}, D195--D202.
#'
#' @source Heinz S, Benner C, Spann N, Bertolino E et al. (2010 May 28) Simple Combinations
#'   of Lineage-Determining Transcription Factors Prime cis-Regulatory
#'   Elements Required for Macrophage and B Cell Identities. Mol Cell, \bold{38(4):576-589}.
#'   PMID: \href{http://www.ncbi.nlm.nih.gov/sites/entrez?Db=Pubmed&term=20513432[UID]}{20513432}
#'
#' @source J Wang, J Zhuang, S Iyer, XY Lin, et al. (2012) Sequence features and
#'   chromatin structure around the genomic regions bound by 119 human transcription
#'   factors. Genome Research, \bold{22 (9)}, 1798-1812, doi:10.1101/gr.139105.112
#'
#' @source Pouya Kheradpour and Manolis Kellis (2013 December 13) Systematic
#'   discovery and characterization of regulatory motifs in ENCODE TF binding
#'   experiments. Nucleic Acids Research, doi:10.1093/nar/gkt1249
#' @examples
#' data(motifbreakR_motif)
#' motifbreakR_motif
"motifbreakR_motif"

#' Example Results from motifbreakR
#'
#' This contains example results from motifbreaker for use in examples from the help docs
#'
#' @format \code{\link[GenomicRanges]{GRanges}} output from \code{motifbreakR}
#'
#' @return \code{\link[GenomicRanges]{GRanges}} object. See \code{\link{motifbreakR}} for information on it's structure.
#'
#' @examples
#' data(example.results)
#' example.results
#'
"example.results"

