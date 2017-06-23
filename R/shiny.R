
#'
#' Main sample overview
#' @noRd
uiMainOverview <- fluidPage(
    column(3,
       h1("FraseR object overview"),
       htmlOutput("fraserShow")
    ),
    column(9,
        h1("Sample Annotation"),
        DT::dataTableOutput("sampleData")
    )
)

#'
#' main result page for shiny
#' @noRd
uiMainResults <- fluidPage(

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
        column(2, textInput("sampleReg", label=h3("Sample search"),
            value="", placeholder="MUC13[46][45] MUC1406")
        ),
        column(2, textInput("symbolReg", label = h3("Gene search"),
            value = "", placeholder = "TIMMDC1 SERA.*")
        )
    ),
    DT::dataTableOutput("results"),
    h1("Links to public resources"),
    DT::dataTableOutput("resultsLinks")

)

#'
#' sample specific graphs and infos
#' @noRd
getUiSampleResults <- function(fds){
    sidebarLayout(
        sidebarPanel = sidebarPanel(width=3,
            selectInput("sampleID", label="Select a sample",
                choices = samples(fds), selected = 1
            )
        ),
        mainPanel = mainPanel(fluidPage(
            plotlyOutput("plotPsi3", height="250px"),
            plotlyOutput("plotPsi5", height="250px"),
            plotlyOutput("plotPsiSite", height="250px"),
            DT::dataTableOutput("selectedPoints")
        ))
    )
}

#'
#' overall main navigation
#' @noRd
getUiNavPage <- function(fds){
    navbarPage("FraseR",
        tabPanel("Main Overview", uiMainOverview),
        tabPanel("All Results", uiMainResults),
        tabPanel("Sample specific", getUiSampleResults(fds)),
        tabPanel("About", p(paste(
            "This shiny app is part of the package 'FraseR'",
            "created by Christian Mertes and Julien Gagneur."))
        )
    )
}

#'
#' main sample results
#' @noRd
getMainPsiTypePanel <- function(sampleID, psiType, fds){
    plotdata <- plotVolcano(fds, sampleID, psiType, source=psiType)
    shinyPlotDF[[psiType]] <<- plotdata[["plotDF"]]
    plotdata[["plot"]] %>% layout(dragmode="select", showlegend=TRUE,
            legend = list(x = 1, y = 0.0, title = "&#936; filter")
    )
}

#'
#' FraseR shiny main server function
#' @noRd
serverMain <- function(input, output) {
    stopifnot(is(shinyFds, "FraseRDataSet"))

    # overview
    output$sampleData <- DT::renderDataTable(
        as.data.table(colData(shinyFds)), escape = FALSE
    )
    output$fraserShow <- renderUI(HTML(
        paste(capture.output(show(shinyFds)), collapse = "</br>")
    ))

    # main results
    filteredResDT <- reactive({
        tmpRes <- shinyFdsRes[
            mcols(shinyFdsRes)$pvalue      <= input$pvalue &
                mcols(shinyFdsRes)$psiValue    <= input$psivalue &
                abs(mcols(shinyFdsRes)$zscore) >= input$zscore
            ]
        try(silent=TRUE, expr={
            tmpRes <- tmpRes[!is.na(mcols(shinyFdsRes)$hgnc_symbol)]
            tmpRes <- tmpRes[
                grepl(perl=TRUE, ignore.case=TRUE,
                      pattern=gsub("\\s+", "|", input$symbolReg, perl=TRUE),
                      x=mcols(tmpRes)$hgnc_symbol
                )
            ]
        })
        tmpRes <- tmpRes[grepl(perl=TRUE, ignore.case=TRUE,
            pattern=gsub("\\s+", "|", input$sampleReg, perl=TRUE),
            x=mcols(tmpRes)$sampleID
        )]
        as.data.table(tmpRes)
    })

    getSelectedLinks <- reactive({
        if(is.null(input$results_rows_selected) |
            length(input$results_rows_selected) < 1){
            return(data.table())
        }

        createFullLinkTable(
            filteredResDT()[input$results_rows_selected],
            TRUE
        )

    })

    output$results <- DT::renderDataTable(filteredResDT())
    if("hgnc_symbol" %in% colnames(mcols(shinyFdsRes))){
        output$resultsLinks <- DT::renderDataTable(
            getSelectedLinks(), escape = FALSE
        )
    }

    # create main figure panel
    output$plotPsi3 <- renderPlotly({
        getMainPsiTypePanel(input$sampleID, "psi3", shinyFds)
    })
    output$plotPsi5 <- renderPlotly({
        getMainPsiTypePanel(input$sampleID, "psi5", shinyFds)
    })
    output$plotPsiSite <- renderPlotly({
        getMainPsiTypePanel(input$sampleID, "psiSite", shinyFds)
    })

    output$selectedPoints <- DT::renderDataTable(escape = FALSE, {

        # Get subset based on selection
        event.data <- lapply(c("psi3", "psi5", "psiSite"), function(x){
            ans <- event_data("plotly_selected", source = x)
            shinyPlotDF[[x]][psiType==x &
                    traceNr %in% ans$curveNumber &
                    pointNr %in% ans$pointNumber
            ]
        })
        # merge data.tables and the get the symbols and create the link table
        pointsOfInterest <- rbindlist(event.data)
        if(nrow(pointsOfInterest) == 0){
            data.table()
        } else {
            pointsOfInterest[,hgnc_symbol:=symbol]
            pointsOfInterest[,seqnames:=chr]
            createFullLinkTable(pointsOfInterest, TRUE)
        }
    })
}


#'
#' Present the FraseR results as shiny app
#'
#' @examples
#'     fds <- createTestFraseRDataSet()
#'
#'     # for interactive sessions
#'     myShinyApp <- FraseRShinyApp(fds)
#'
#'     # for running a shiny application as server
#'     myShinyAppObj <- FraseRShinyApp(fds, server=TRUE)
#'
#' @export
FraseRShinyApp <- function(fds, fdsres=NULL, server=!interactive(), ...){
    if(is.null(fdsres)){
        fdsres <- results(fds)
    }

    shinyFds    <<- fds
    shinyFdsRes <<- fdsres
    shinyPlotDF <<- list(psi3=NULL, psi5=NULL, psiSite=NULL)

    shinyObj <- list(
        shinyFds    = shinyFds,
        shinyFdsRes = shinyFdsRes,
        serverMain  = serverMain,
        uiMain      = getUiNavPage(shinyFds)
    )

    if(!server){
       shinyApp(
            ui = shinyObj$uiMain,
            server = shinyObj$serverMain,
            ...
        )
    } else {
        return(list(
            shinyObj=shinyObj,
            shinyFds=shinyFds,
            shinyFdsRes=shinyFdsRes
        ))
    }
}
