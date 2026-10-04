# library(BiocManager)
# options(repos = BiocManager::repositories())
# library(shiny)
# library(bslib)
# library(bsicons)
# library(BSgenome)
# library(motifbreakR)
# library(MotifDb)
# library(biomaRt)
# library("BSgenome.Hsapiens.1000genomes.hs37d5")
# library("BSgenome.Hsapiens.UCSC.hg19")
# library("BSgenome.Hsapiens.UCSC.hg19.masked")
# library("BSgenome.Hsapiens.NCBI.GRCh38")
# library("BSgenome.Hsapiens.UCSC.hg38.masked")
# library("BSgenome.Hsapiens.UCSC.hg38")

strToObj <- function(x) {eval(parse(text=x))}

ui <- bslib::page_navbar(window_title = "motifbreakR", fillable = FALSE,
                         theme = bslib::bs_theme(bootswatch = "sketchy"),
                         title = "motifbreakR",
                         id = "main",
                         sidebar = bslib::sidebar(title = "Input Data",
                                                  width = 450,
                                                  bslib::card(
                                                    selectInput(
                                                      inputId = "selected.genome",
                                                      label = "Select Genome",
                                                      choices = c("Choose BSgenome" = "",
                                                                  c("BSgenome.Hsapiens.1000genomes.hs37d5",
                                                                    "BSgenome.Hsapiens.UCSC.hg19",
                                                                    "BSgenome.Hsapiens.UCSC.hg19.masked",
                                                                    "BSgenome.Hsapiens.NCBI.GRCh38",
                                                                    "BSgenome.Hsapiens.UCSC.hg38.masked",
                                                                    "BSgenome.Hsapiens.UCSC.hg38"))
                                                    ),
                                                    splitLayout(
                                                      radioButtons(inputId = "selected.format",
                                                                   label = "Format",
                                                                   choices = c("rsID", "BED", "VCF"),
                                                                   selected = "rsID"),
                                                      uiOutput("file_input_data"),
                                                      cellWidths = c("25%", "75%")
                                                    ),
                                                    actionButton(
                                                      inputId = "submit.analysis",
                                                      label = "Submit"
                                                    ),
                                                  ),
                                                  bslib::card(
                                                    selectInput(
                                                      inputId = "motif.set",
                                                      label = "Select Motif List",
                                                      choices = c("Choose Motif Set" = "", unique(mcols(MotifDb)$dataSource)),
                                                      multiple = TRUE,
                                                      selectize = TRUE
                                                    ),
                                                    numericInput(
                                                      inputId = "selected.pval",
                                                      label = "P-value threshold",
                                                      value = 5e-5,
                                                      max = 1, min = 0, step = 1e-6,
                                                    ),
                                                    bslib::tooltip(
                                                      span(numericInput(
                                                        inputId = "gran.pval",
                                                        label = list("P-value granularity", bsicons::bs_icon("question-circle-fill")),
                                                        value = 1e-4,
                                                        max = 1, min = 0, step = 1e-7,
                                                      )),
                                                      "The granularity of the the PWMS (basically rounding)\nused for calculating the p-value.\n1e-4 for speed, 1e-7 for greater accuracy",
                                                      placement = "right"
                                                    ),
                                                    actionButton(inputId = "execute.analysis", label = "Run")
                                                  )),
                         bslib::nav_panel(title = "Input",
                                   fillable = FALSE,
                                   bslib::card(
                                     bslib::card_header(
                                       class = "bg-dark",
                                       "Input SNVs"),
                                     DT::DTOutput("user.input")),
                                   bslib::card(
                                     bslib::card_header(
                                       class = "bg-dark",
                                       "Individual Motifs"),
                                     DT::DTOutput("motif.selections"))
                         ),
                         bslib::nav_panel(title = "Results",
                                   splitLayout(p("Results may be plotted for a table filtered to include a single variant id, with a single alternate allele"),
                                               uiOutput("downloadData"),
                                               cellArgs = list(style='white-space: normal;')),
                                   shiny::selectInput(inputId = "select_results",
                                                      label = "Selected columns to display", "",
                                                      width = '100%',
                                                      multiple = TRUE),
                                   DT::DTOutput("motifbreakr.results"),
                                   actionButton(inputId = "plot.motifbreakR", label = "Plot Displayed Results", disabled = TRUE, width = "100%")
                         ),
                         bslib::nav_panel(title = "Plot",
                                   splitLayout(textOutput("text"),downloadButton("downloadPlot", "Download Plot")),
                                   splitLayout(numericInput("plotwidth", "Width:", min = 200, max = 3000, value = 750, step = 5),
                                               numericInput("plotheight", "Height:", min = 200, max = 3000, value = 500, step = 5)),
                                   plotOutput("motifbreakr.plot")),
                         bslib::nav_spacer(),
                         bslib::nav_item(list(bslib::input_dark_mode(), " set dark mode"), align= "right")

)

server <- function(input, output, session) {

  output$file_input_data <- renderUI({
    req(input$selected.format)
    switch(input$selected.format,
           rsID = textAreaInput(inputId = "user.data",
                                label = "Enter rsid values, one on each line"),
           BED = fileInput(inputId = "user.data",
                           label = "Upload SNPs",
                           accept = ".bed"),
           VCF = fileInput(inputId = "user.data",
                           label = "Upload SNPs",
                           accept = c(".vcf", ".vcf.gz"))
    )
  })

  user.input <- eventReactive(input$submit.analysis, {
    req(input$user.data, input$selected.genome)
    withProgress(message = "Data getting ready",
                 detail = "This can take a while, looking up variants ...",
                 value = 0, {
                   inFile <- input$user.data
                   if (is.null(inFile)) {
                     return(NULL)
                   }

                   if(input$selected.format == "rsID") {
                     inFile <- strsplit(inFile, "\n")[[1]]
                   }

                   incProgress(1/3)

                   userBSgenome <- strToObj(paste(input$selected.genome, input$selected.genome, sep="::"))

                   if(organism(userBSgenome) == "Homo sapiens" & (input$selected.format %in% c("rsID", "BED"))) {
                     genome_version <- switch (input$selected.genome,
                                               BSgenome.Hsapiens.1000genomes.hs37d5 = "GRCh37",
                                               BSgenome.Hsapiens.UCSC.hg19 = "GRCh37",
                                               BSgenome.Hsapiens.UCSC.hg19.masked = "GRCh37",
                                               BSgenome.Hsapiens.NCBI.GRCh38 = NULL,
                                               BSgenome.Hsapiens.UCSC.hg38.masked = NULL,
                                               BSgenome.Hsapiens.UCSC.hg38 = NULL,
                                               "no_rsid")
                     if(isTRUE(genome_version == "no_rsid")) {
                       renderPrint(stop("rsID values are not accepted for this genome version\n Please select hg19 or hg38 compatible genomes."))
                     }
                     bmsnp <- motifbreakR:::loadMartObj(genome_version)
                   } else {
                     if(input$selected.format == "rsID") renderPrint(stop("motifbreakR doesn't accept rsid values for non human genomes"))
                   }
                   incProgress(1/3)

                   user.input <- switch(input$selected.format,
                                        rsID = snps.from.rsid(inFile,
                                                              biomart.dataset = bmsnp,
                                                              search.genome = userBSgenome),
                                        BED = snps.from.file(file = inFile$datapath, search.genome = userBSgenome,
                                                             biomart.dataset = bmsnp, format = "bed"),
                                        VCF = snps.from.file(file = inFile$datapath, search.genome = userBSgenome,
                                                             format = "vcf"))

                   incProgress(1/3)
                   return(user.input)
                 })
  })

  mList <- reactive({
    mList <- input$motif.set
    if (is.null(mList)) {
      return(NULL)
    }
    mList <- MotifDb[mcols(MotifDb)$dataSource %in% input$motif.set]
    return(mList)
  })

  motifFilter <- reactiveValues(rowind = NULL)

  observeEvent(input$motif.selections_rows_all, {
    motif_displayedRowsIndices <- input$motif.selections_rows_all
    motifFilter$rowind <- motif_displayedRowsIndices
  })

  observeEvent(input$execute.analysis, {
    bslib::nav_select("main", "Results")
  })


  results <- eventReactive(input$execute.analysis, {
    withProgress(message = "Running Analysis:\n",
                 detail = "running motifbreakR",
                 value = 0, {
                   incProgress(1/3)
                   res <- motifbreakR(snpList = user.input(),
                                      pwmList = mList()[motifFilter$rowind],
                                      threshold = input$selected.pval,
                                      method = "ic",
                                      filterp = TRUE)
                   incProgress(1/3, detail = "Finding Matching Peaks")
                   genome_version <- genome(user.input())[[1]]
                   res <- findSupportingRemapPeaks(res,
                                                   genome = genome_version,
                                                   TFClass = TRUE)
                   incProgress(1/3, detail = "Calculating Pvalue")
                   res <- calculatePvalue(res, granularity = input$gran.pval)
                   return(res)
                   }
                 )})

  output$user.input <- DT::renderDT({
    req(user.input())
    input.df <- as.data.frame(user.input(), row.names = NULL)
    input.df <- input.df[, c("seqnames", "start", "end", "SNP_id", "REF", "ALT")]
    DT::datatable(input.df, filter = "top", rownames = FALSE)
  })

  results.colnames <- c("seqnames",
                        "start",
                        "end",
                        "width",
                        "strand",
                        "SNP_id",
                        "REF",
                        "ALT",
                        "varType",
                        "windowIdxRef",
                        "motifStrandRef",
                        "windowIdxAlt",
                        "motifStrandAlt",
                        "geneSymbol",
                        "dataSource",
                        "providerName",
                        "providerId",
                        "seqMatch",
                        "pwmConsensus",
                        "pctRef",
                        "pctAlt",
                        "scoreRef",
                        "scoreAlt",
                        "pValueRef",
                        "pValueAlt",
                        "strongerIn",
                        "alleleDiff",
                        "alleleEffectSize",
                        "effect",
                        "matchingBindingEvent",
                        "pValueEffect")

  output$motifbreakr.results <- DT::renderDT({
    empty <- as.data.frame(matrix(ncol = 6,
                                  nrow = 1,
                                  dimnames = list(c(),
                                                  c("seqnames",
                                                    "start",
                                                    "end",
                                                    "width",
                                                    "strand",
                                                    "SNP_id",
                                                    "REF",
                                                    "ALT",
                                                    "varType",
                                                    "windowIdxRef",
                                                    "motifStrandRef",
                                                    "windowIdxAlt",
                                                    "motifStrandAlt",
                                                    "geneSymbol",
                                                    "dataSource",
                                                    "providerName",
                                                    "providerId",
                                                    "seqMatch",
                                                    "pwmConsensus",
                                                    "pctRef",
                                                    "pctAlt",
                                                    "scoreRef",
                                                    "scoreAlt",
                                                    "pValueRef",
                                                    "pValueAlt",
                                                    "strongerIn",
                                                    "alleleDiff",
                                                    "alleleEffectSize",
                                                    "effect",
                                                    "matchingBindingEvent",
                                                    "pValueEffect"))))
    DT::datatable(empty, filter = "top", rownames = F)
    })

  observe({
    req(user.input())
    inputdata.colnames <- names(as.data.frame(user.input(), row.names = NULL))
    updateSelectInput(session, "userInputNames",
                      label = NULL,
                      choices = inputdata.colnames,
                      selected = c("seqnames", "start", "end", "SNP_id", "REF", "ALT"))
  })

  output$motif.selections <- DT::renderDT({
    req(mList())
    motif.df <- as.data.frame(mcols(mList()))
    motif.df$dataSource <- as.factor(motif.df$dataSource)
    motif.df$geneSymbol <- as.factor(motif.df$geneSymbol)
    motif.df$organism <- factor(motif.df$organism)
    motif.df$tfFamily <- as.factor(motif.df$tfFamily)
    motif.df$experimentType <- as.factor(motif.df$experimentType)
    motif.df$bindingDomain <- as.factor(motif.df$bindingDomain)
    motif.df$bindingSequence <- as.factor(motif.df$bindingSequence)
    DT::datatable(motif.df,
                  filter = 'top',
                  rownames = FALSE)
  })

  updateSelectInput(session, "select_results",
                    label = NULL,
                    choices = results.colnames,
                    selected = c("seqnames", "start", "end", "SNP_id", "REF", "ALT",
                                 "geneSymbol", "pctRef", "pctAlt", "scoreRef",
                                 "scoreAlt", "alleleDiff", "effect"))


  output$motifbreakr.results <- DT::renderDT({
    req(results())
    subset_data <- results()
    subset_data <- as.data.frame(subset_data, row.names = NULL)
    subset_data$SNP_id <- as.factor(subset_data$SNP_id)
    subset_data$geneSymbol <- as.factor(subset_data$geneSymbol)
    subset_data$geneSymbol <- as.factor(subset_data$geneSymbol)
    subset_data$providerName <- as.factor(subset_data$providerName)
    subset_data$providerId <- as.factor(subset_data$providerId)
    subset_data$effect <- as.factor(subset_data$effect)
    subset_data$pValueEffect <- as.factor(subset_data$pValueEffect)
    subset_data <- subset_data[, input$select_results, drop = FALSE]
    subset_data_cols <- c("pctRef", "pctAlt", "scoreRef", "scoreAlt", "alleleDiff")
    subset_data_cols <- subset_data_cols[subset_data_cols %in% input$select_results]
    if("matchingBindingEvent" %in% c(input$select_results)) {
      convert_to_char <- which(is.na(subset_data[, "matchingBindingEvent"]), arr.ind = TRUE)
      subset_data <- DataFrame(subset_data)
      match_peak_col_i_o <- which(names(subset_data) == c("matchingBindingEvent"))
      no_na <- lapply(subset_data[, match_peak_col_i_o],
                      function(x) {
                        ifelse(is.na(x), x <- "", x)
                      }
      )
      subset_data[["matchingBindingEvent"]] <- no_na
    }
    subset_data <- as.data.frame(subset_data)
    DT::formatRound(DT::datatable(subset_data, selection = "multiple", filter = "top", rownames = F,
                                  options = list(
                                    rowCallback = DT::JS(
                                      "function(row, data) {",
                                      "for (i = 1; i < data.length; i++) {",
                                      "if (data[i] === '' && Math.abs(data[i])<0.01){",
                                      "$('td:eq('+i+')', row).html(data[i].toExponential(1));",
                                      "}",
                                      "}",
                                      "}"))
                                  ),
      columns = subset_data_cols, digits = 3)
  })

  output$downloadData <- renderUI({
    req(results())
    downloadButton("downloadableData", "Download Complete Data")
  })

  output$downloadableData <- downloadHandler(
    filename = function() { "mb_output.csv" },
    content = function(file) {
      exportMBtable(results = results(), file = file, format = "csv")
    }
  )

  plotFilter <- reactiveValues(rowind = NULL, snp = NULL, alt = NULL)

  observeEvent(input$motifbreakr.results_rows_all, {
    displayedRowsIndices <- input$motifbreakr.results_rows_all
    filteredData <- results()[displayedRowsIndices, ]
    snp_to_plot <- NULL
    snp_alt <- NULL
    if(length(unique(filteredData$SNP_id)) == 1 & length(unique(filteredData$ALT))) {
      snp_to_plot <- as.character(filteredData$SNP_id[[1]])
      snp_alt <- as.character(filteredData$ALT[[1]])
      updateActionButton(session, "plot.motifbreakR",
                         disabled = FALSE)
    } else {
      updateActionButton(session, "plot.motifbreakR",
                         disabled = TRUE)
    }
    plotFilter$rowind <- displayedRowsIndices
    plotFilter$snp <- snp_to_plot
    plotFilter$alt <- snp_alt
  })

  observeEvent(input$plot.motifbreakR, {
    bslib::nav_select("main", "Plot")
    output$text <- renderText(paste0("Rendering Plot for ", plotFilter$snp, " ",
                                     "with alt allele ", plotFilter$alt, "."))
    output$motifbreakr.plot <- renderPlot({
      withProgress(message = "Plotting Result: ",
                   detail = "selecting variant",
                   value = 1/3, {
      plot_data <- results()[plotFilter$rowind]
      snp_to_plot <- plotFilter$snp
      incProgress(1/3, detail = "plotting variant")
      plotMB(plot_data, rsid = snp_to_plot)
      incProgress(1/3)})
      }, height = function() input$plotheight, width = function() input$plotwidth)
    })

  output$downloadPlot <- downloadHandler(
    filename = function() { "mb_output.pdf" },
    content = function(file) {
      plot_data <- results()[plotFilter$rowind]
      snp_to_plot <- plotFilter$snp
      pdf(file,
          width = input$plotwidth/100,
          height = input$plotheight/100)
      plotMB(plot_data, rsid = snp_to_plot)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
