library(shiny)
library(BSgenome)
library(motifbreakR)
library(MotifDb)
library(shinythemes)
#data(encodemotif); data(factorbook)

#MotifDb <- c(MotifDb, encodemotif, factorbook)

strToObj <- function(x) {eval(parse(text=x))}

ui <- fluidPage(
    titlePanel("MotifBreakR"),
    sidebarLayout(
        sidebarPanel(
            wellPanel(
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
                    max = 1, min = 0, step = 1e-6)
            )
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Input Data",
                    navbarPage("User Input",
                        tabPanel("Input SNVs", DT::dataTableOutput("user.input")),
                        tabPanel("Individual Motifs", DT::dataTableOutput("motif.selections"))
                    )
                ),
                tabPanel("Results", DT::dataTableOutput("motifbreakr.results")),
                tabPanel("Plot", "This is where the plots will be")
            )
        )
    )
)

server <- function(input, output) {
    user.input <- reactive({
        withProgress(message = "Data getting ready",
                     detail = "This can take a while, looking up variants ...",
                     value = 0, {
                         inFile <- input$user.data

                         if (is.null(inFile)) {
                             return(NULL)
                         }

                         inFile <- readr::read_tsv(file = inFile$datapath, col_names = FALSE)
                         incProgress(1/3)
                         if (!requireNamespace(input$selected.snplocs, quietly = TRUE, character.only = TRUE)) {
                             stop("Cannot locate SNPlocs package ", input$selected.snplocs)
                         } else {
                             userSNPlocs <- eval(parse(text = paste(input$selected.snplocs, input$selected.snplocs, sep="::")))
                         }

                         if(!requireNamespace(input$selected.genome, quietly = TRUE, character.only = TRUE)) {
                             stop("Cannot locate BSgenome package ", input$selected.genome)
                         } else {
                             userBSgenome <- eval(parse(text = paste(input$selected.genome, input$selected.genome, sep="::")))
                         }
                         incProgress(1/3)
                         user.input <- switch(input$selected.format,
                                              rsID = snps.from.rsid(inFile$X1,
                                                                    dbSNP = userSNPlocs,
                                                                    search.genome = userBSgenome),
                                              BED = snps.from.file(inFile, format = "blah"),
                                              VCF = snps.from.file(inFile, format = "blah"))
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

    output$user.input <- DT::renderDataTable({
        # if(is(user.input(), "Observer")) return(NULL)
        # message(class(user.input))
        # message("is null? ", is.null(user.input))
        my_ui <- observe({as.data.frame(user.input())})
        as.data.frame(my_ui)
    })

    output$motif.selections <- DT::renderDataTable({
    # output$motif.selections <- renderPrint({
        # message(class(mList))
        # message("is null? ", is.null(mList))
        as.data.frame(mcols(observe({mList})))
    })

    output$motifbreakr.results <- DT::renderDataTable({

        if (!exists(mList) | is.null(mList)) {
            return(NULL)
        }

        results <- motifbreakR(snpList = user.input,
                               pwmList = mList,
                               threshold = input$selected.pval,
                               method = "ic",
                               filterp = TRUE)
        as.data.frame(results)
    })
}

shinyApp(ui = ui, server = server)

