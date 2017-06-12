
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
        sidebarPanel = sidebarPanel(width=4,
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
#' FraseR shiny main server function
#' @noRd
serverMain <- function(input, output) {
    stopifnot(is(shinyFds, "FraseRDataSet"))

    # overview
    output$sampleData <- DT::renderDataTable(as.data.table(colData(shinyFds)))
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

    # main sample results
    output$plotPsi3 <- renderPlotly({
        plotdata <- plotVolcano(fds, input$sampleID, "psi3", source="psi3")
        shinyPlotted_psi3 <<- plotdata[["pointPlotted_psi3"]]
        plotdata[["plot"]] %>% layout(dragmode="select", showlegend=TRUE,
                legend = list(x = 1, y = 0.2, title = "&#936; filter")
        )
    })
    output$plotPsi5 <- renderPlotly({
        plotdata <- plotVolcano(fds, input$sampleID, "psi5", source="psi5")
        shinyPlotted_psi5 <<- plotdata[["pointPlotted_psi5"]]
        plotdata[["plot"]] %>% layout(dragmode="select", showlegend=TRUE,
                legend = list(x = 1, y = 0.2, title = "&#936; filter")
        )
    })
    output$plotPsiSite <- renderPlotly({
        plotdata <- plotVolcano(fds, input$sampleID, "psiSite", source="psiSite")
        shinyPlotted_psiSite <<- plotdata[["pointPlotted_psiSite"]]
        plotdata[["plot"]] %>% layout(dragmode="select", showlegend=TRUE,
                legend = list(x = 1, y = 0.2, title = "&#936; filter")
        )
    })

    output$selectedPoints <- DT::renderDataTable({

        # Get subset based on selection
        event.data <- lapply(c("psi3", "psi5", "psiSite"), function(x){
            ans <- event_data("plotly_selected", source = x)
            as.data.table(ans)
        })
        event.data <- rbindlist(event.data)

        data.table(event.data)
        #data.table(subset(
        #       event.data, curveNumber == 0)$pointNumber + 1
        #)

#
#         # Get number of malignant and benign cases from selection
#         malig.class <- subset(plot.df, Class == "malignant")[subset(event.data, curveNumber == 0)$pointNumber + 1,]
#         benign.class <- subset(plot.df, Class == "benign")[subset(event.data, curveNumber == 1)$pointNumber + 1,]
#
#         # Combine
#         plot.subset <- rbind(malig.class, benign.class)
#
#         # Summarize
#         plot.summ <- plot.subset %>%
#             group_by(x, y, Class) %>%
#             summarize(Count = n())
#
#         # Assign to parent frame
#         plot.summ <<- plot.summ
#
#         # Plot
#         plot_ly(plot.summ, x = ~Class, y = ~Count, type = "bar", source = "select", color = ~Class) %>%
#             layout(title = "No. of Malignant and Benign cases <br> in Selection",
#                    plot_bgcolor = "6A446F",
#                    yaxis = list(domain = c(0, 0.9)))
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
