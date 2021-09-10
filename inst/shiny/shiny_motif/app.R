library(shiny)
library(BSgenome)
library(motifbreakR)
library(MotifDb)
library(shinythemes)
#data(encodemotif); data(factorbook)

#MotifDb <- c(MotifDb, encodemotif, factorbook)

strToObj <- function(x) {eval(parse(text=x))}

ui <- fluidPage(theme = shinytheme("simplex"),
    titlePanel("MotifBreakR"),
    sidebarLayout(
        sidebarPanel(
            wellPanel(
                radioButtons(inputId = "indels",
                             label = "Variant Types:",
                             choices = c("SNVs and/or Indels", "SNVs only")),
                selectInput(inputId = "selected.genome",
                    label = "Select Genome",
                    choices = c("Choose BSgenome" = "", available.genomes())),
                selectInput(inputId = "selected.snplocs",
                    label  = "Select SNP Catalogue",
                    choices = c("Choose SNPlocs" = "", available.SNPs())),
                splitLayout(
                    radioButtons(inputId = "selected.format",
                        label = "Format",
                        choices = c("rsID", "BED", "VCF"),
                        selected = "rsID"),
                    fileInput(inputId = "user.data",
                        label = "Upload SNPs"),
                    cellWidths = c("25%", "75%")
                )
            ),
            wellPanel(
                selectInput(inputId = "motif.set",
                    label = "Select Motif List",
                    choices = c("Choose Motif Set" = "", unique(mcols(MotifDb)$dataSource)),
                    multiple = TRUE,
                    selectize = TRUE),
                numericInput(inputId = "selected.pval",
                    label = "P-value threshold",
                    value = 5e-5,
                    max = 1, min = 0, step = 1e-6),
                actionButton(inputId = "execute.analysis", label = "Run")
            ),
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Input Data",
                    navbarPage("User Input",
                        tabPanel("Input SNVs", DT::dataTableOutput("user.input")),
                        tabPanel("Individual Motifs", DT::dataTableOutput("motif.selections"))
                    )
                ),
                tabPanel("Results",
                         selectInput("select_results", "Select columns to display",
                                     choices = c("seqnames", "start", "end", "width", "strand",
                                       "SNP_id", "REF", "ALT", "varType", "motifPos",
                                       "geneSymbol", "dataSource", "providerName", "providerID",
                                       "seqMatch", "pctRef", "pctAlt", "scoreRef", "scoreAlt",
                                       "Refpvalue", "Altpvalue", "altPos", "alleleDiff", "effect"),
                                     selected = c("seqnames", "start", "end", "SNP_id", "REF", "ALT",
                                                  "geneSymbol", "dataSource", "providerName",
                                                  "pctRef", "pctAlt", "scoreRef", "scoreAlt", "alleleDiff", "effect"),
                                     width = '80%',
                                     multiple = TRUE),
                         DT::dataTableOutput("motifbreakr.results")),
                tabPanel("Plot", "This is where the plots will be")
            )
        )
    )
)

server <- function(input, output) {
    output$dbsnp <- renderUI({
        if (input$indels == "SNVs only") {

        }
    })
    user.input <- reactive({
        req(input$user.data, input$selected.snplocs, input$selected.genome)
        withProgress(message = "Data getting ready",
                     detail = "This can take a while, looking up variants ...",
                     value = 0, {
                         inFile <- input$user.data

                         if (is.null(inFile)) {
                             return(NULL)
                         }

                         inFiler <- readr::read_tsv(file = inFile$datapath, col_names = FALSE)
                         incProgress(1/3)

                         if (input$selected.snplocs == "") {
                             userSNPlocs <- NULL
                         } else {
                             userSNPlocs <- eval(parse(text = paste(input$selected.snplocs, input$selected.snplocs, sep="::")))
                         }
                         userBSgenome <- eval(parse(text = paste(input$selected.genome, input$selected.genome, sep="::")))

                         incProgress(1/3)
                         user.input <- switch(input$selected.format,
                                              rsID = snps.from.rsid(inFiler$X1,
                                                                    dbSNP = userSNPlocs,
                                                                    search.genome = userBSgenome),
                                              BED = snps.from.file(file = inFile$datapath, search.genome = userBSgenome,
                                                                   dbSNP = userSNPlocs, format = "bed", indels = ifelse(input$indels == "SNVs only", FALSE, TRUE)),
                                              VCF = snps.from.file(file = inFile$datapath, search.genome = userBSgenome,
                                                                   dbSNP = userSNPlocs, format = "vcf", indels = ifelse(input$indels == "SNVs only", FALSE, TRUE)))
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

    results <- eventReactive(input$execute.analysis, {
        motifbreakR(snpList = user.input(),
                    pwmList = mList(),
                    threshold = input$selected.pval,
                    method = "ic",
                    legacy.score = FALSE,
                    filterp = TRUE)
    })

    output$user.input <- DT::renderDataTable({
        req(user.input())
        as.data.frame(user.input())
    })

    output$motif.selections <- DT::renderDataTable({
        as.data.frame(mcols(mList()))
    })


    output$motifbreakr.results <- DT::renderDataTable({
        DT::formatRound(DT::datatable(as.data.frame(results(), row.names = NULL)[, input$select_results]),
                                      columns = c("pctRef", "pctAlt", "scoreRef", "scoreAlt", "alleleDiff"), digits = 3)
    })
}

shinyApp(ui = ui, server = server)

