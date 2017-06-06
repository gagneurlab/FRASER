#'
#' creates the shiny app object
#'
#' @noRd
FraseRshinyApp <- function(shinyFds, shinyFdsRes){
    require(FraseR)

#'
#' Main sample overview
#'
uiMainOverview <- fluidPage(
    titlePanel("Main Overview"),

    h1("The FraseR result shiny server"),

    DT::dataTableOutput("sampleData")
)

#'
#' main result page for shiny
#' @noRd
uiMainResults <- fluidPage(
    titlePanel("FraseR main results"),

    fluidRow(
        column(2, numericInput("pvalue", label=h3("Max P-value"),
            value=1e-5, min=0, max=1)
        ),
        column(2, numericInput("zscore", label=h3("Min Z-score"),
            value=2.5, min=0, max=50)
        ),
        column(2, numericInput("psivalue", label=h3("Max PSI-value"),
            value=0.3, min=0, max=1)
        ),
        column(3, textInput("sampleReg", label=h3("Sample search"),
            value="", placeholder="MUC13[46][45] MUC1406")
        ),
        column(3, textInput("symbolReg", label = h3("Gene search"),
            value = "", placeholder = "TIMMDC1 SERA.*")
        )
    ),
    DT::dataTableOutput("results")
)

#'
#'
#'
uiSampleResults <- fluidPage(
    titlePanel("FraseR sample results"),

    sidebarLayout(
        sidebarPanel = sidebarPanel(
            column(3, selectInput("sampleID", label="Select a sample",
                    choices = FraseR::samples(shinyFds), selected = 1
            ))
        ),
        mainPanel = mainPanel(plotlyOutput("plotMainFraseR", height="100%", width="100%"))
    )
)

#'
#'
#'
navPage <- navbarPage("FraseR",
    tabPanel("Main Overview", uiMainOverview),
    tabPanel("All Results", uiMainResults),
    tabPanel("Sample specific", uiSampleResults),
    tabPanel("About", p(paste("This shiny app is part of the package 'FraseR'",
            "created by Christian Mertes and Julien Gagneur."))
    )
)

#'
#' FraseR shiny main server function
#' @noRd
serverMain <- function(input, output) {
    # main overview
    output$sampleData <- DT::renderDataTable(as.data.table(
        colData(shinyFds)
    ))

    # main results
    output$results <- DT::renderDataTable(as.data.table({
        tmpRes <- shinyFdsRes[
            mcols(shinyFdsRes)$pvalue      <= input$pvalue &
            mcols(shinyFdsRes)$psiValue    <= input$psivalue &
            abs(mcols(shinyFdsRes)$zscore) >= input$zscore &
            !is.na(mcols(shinyFdsRes)$hgnc_symbol)
        ]
        tmpRes <- tmpRes[grepl(perl=TRUE, ignore.case=TRUE,
           pattern=gsub("\\s+", "|", input$symbolReg, perl=TRUE),
           x=mcols(tmpRes)$hgnc_symbol
        )]
        tmpRes <- tmpRes[grepl(perl=TRUE, ignore.case=TRUE,
           pattern=gsub("\\s+", "|", input$sampleReg, perl=TRUE),
           x=mcols(tmpRes)$sampleID
        )]
        tmpRes
    }))

    # main sample results
    output$plotMainFraseR <- renderPlotly(
            FraseR:::createMainPlotFraseR(shinyFds, input$sampleID)
    )
}

    return(list(
        shinyFds = shinyFds,
        shinyFdsRes = shinyFdsRes,
        serverMain = serverMain,
        navPage = navPage,
        uiSampleResults = uiSampleResults,
        uiMainResults =uiMainResults
    ))
}

#'
#' Present the FraseR results in a shiny app
#'
#' @export
FraseRShiny <- function(fds, fdsres=NULL, options=list(), server=!interactive()){
    if(is.null(fdsres)){
        shinyFdsRes <- FraseR::results(fds)
    } else {
        shinyFdsRes <- fdsres
    }
    shinyFds <- fds

    shinyObj <- FraseRshinyApp(fds, shinyFdsRes)

    if(!server){
        shinyApp(
            ui = shinyObj$navPage,
            server = shinyObj$serverMain,
            options = options
        )
    } else {
        return(list(
            shinyObj=shinyObj,
            shinyFds=shinyFds,
            shinyFdsRes=shinyFdsRes
        ))
    }
}
