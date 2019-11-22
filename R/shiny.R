
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
            ),
            plotOutput("plotPsiDist"),
            plotOutput("plotRawCountDist"),
            plotOutput("plotZscoreDist")
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
        tabPanel("About", uiAbout())
    )
}

uiAbout <- function(){
    fluidPage(
        fluidRow(
            h3(paste(
                "This shiny app is part of the package 'FraseR'",
                "created by Christian Mertes and Julien Gagneur."
            ))
        ),
        fluidRow(
            h3("Session Info"),
            p(paste("This FraseR was compiled on: ",
                    as.Date.POSIXct(min(unlist(
                        file.info(system.file(package = "FraseR"))[
                            paste0(c('c', 'm', 'a'),'time')
                        ])
                    ))
            ))
        ),
        fluidRow(
            htmlOutput("FraseRSessionInfo")
        )
    )
}

#'
#' main sample results
#' @noRd
getMainPsiTypePanel <- function(sampleID, psiType, fds){
    plotdata <- plotVolcano(fds, sampleID, psiType)
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
        # browser()

        # convert to data.table
        tmpRes <- as.data.table(shinyFdsRes)

        # filter default values
        tmpRes <- tmpRes[
                pvalue      <= input$pvalue   &
                psiValue    <= input$psivalue &
                abs(zscore) >= input$zscore]

        # filter gene names
        if("hgnc_symbol" %in% colnames(tmpRes)){
            tmpRes <- tmpRes[!is.na(hgnc_symbol)]
            tmpRes <- tmpRes[grepl(
                    x=hgnc_symbol, perl=TRUE, ignore.case=TRUE,
                    pattern=gsub("\\s+", "|", input$symbolReg, perl=TRUE))]
        }

        # filter samples
        tmpRes <- tmpRes[grepl(
                x=sampleID, perl=TRUE, ignore.case=TRUE,
                pattern=gsub("\\s+", "|", input$sampleReg, perl=TRUE))]

    })

    getSelectedLinks <- reactive({
        if(is.null(input$results_rows_selected) |
            length(input$results_rows_selected) < 1){
            return(data.table())
        }

        createFullLinkTable(filteredResDT()[input$results_rows_selected], TRUE)

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

    # Get subset based on selection
    selectedPointsAsDT <- reactive({
        event.data <- lapply(c("psi3", "psi5", "psiSite"), function(x){
            ans <- event_data("plotly_selected", source = x)
            shinyPlotDF[[x]][psiType==x &
                    traceNr %in% ans$curveNumber &
                    pointNr %in% ans$pointNumber
            ]
        })

        # merge data.tables and the get the symbols and create the link table
        if(is.null(event.data) || length(event.data) == 0){
            return(data.table())
        }
        return(rbindlist(event.data))
    })

    output$selectedPoints <- DT::renderDataTable(escape = FALSE, {
        #browser()
        pointsOfInterest <- selectedPointsAsDT()

        if(nrow(pointsOfInterest) == 0){
            data.table()
        } else {
            pointsOfInterest[,hgnc_symbol:=symbol]
            pointsOfInterest[,seqnames:=chr]
            createFullLinkTable(pointsOfInterest, TRUE)
        }
    })

    output$plotPsiDist <- renderPlot({
        getForSelectedPointsDistPlot(input$selectedPoints_rows_selected,
                selectedPointsAsDT(), shinyFds,
                dist="psi", sampleID=input$sampleID
        )
    })

    output$plotRawCountDist <- renderPlot({
        getForSelectedPointsDistPlot(input$selectedPoints_rows_selected,
                selectedPointsAsDT(), shinyFds,
                dist="counts", sampleID=input$sampleID
        )
    })

    output$plotZscoreDist <- renderPlot({
        getForSelectedPointsDistPlot(input$selectedPoints_rows_selected,
                selectedPointsAsDT(), shinyFds,
                dist="ocounts", sampleID=input$sampleID
        )
    })

    output$FraseRSessionInfo <- renderPrint({
        capture.output(sessionInfo())
    })
}

getForSelectedPointsDistPlot <- function(selected, selectedFromDT, fds, ...){
    if(is.null(selected)){
        return(NULL)
    }
    ans <- getIdxForSelectedPoints(selected, selectedFromDT, fds)
    if(length(ans) != 2){
        return(NULL)
    }
    getDisttributionPlot(fds, idx = ans[["idx"]],
            psiType = ans[["psiType"]], ...
    )
}

getIdxForSelectedPoints <- function(selected, selectedFromDT, fds){
    if(is.null(selected)){
        return(logical(0))
    }

    selectedData <- selectedFromDT[selected[1]]
    psiType <- selectedData[,psiType]
    querygr <- makeGRangesFromDataFrame(selectedData[,.(chr,start,end)])
    if(psiType=="psiSite"){
        subjgr <- granges(nonSplicedReads(fds))
    } else {
        subjgr <- granges(fds)
    }

    idx <- to(findOverlaps(querygr, subjgr, type="equal"))

    return(list(idx=idx, psiType=psiType))
}

getDisttributionPlot <- function(fds, psiType, dist, idx, sampleID=NULL, ...){

    col <- rep("gray", dim(fds)[2])
    dist <- switch(dist,
            psi=psiType,
            counts=ifelse(psiType=="psiSite", "rawCountsSS", "rawCountsJ"),
            ocounts=paste0("rawOtherCounts_", psiType),
            zscore=paste0("zscore_", psiType)
    )
    data <- as.vector(assays(fds)[[dist]][idx,])
    names(data) <- samples(fds)
    data <- sort(data)
    if(!is.null(sampleID) && any(sampleID %in% names(data))){
        message("color sample: ", sampleID, " on idx: ", which(sampleID %in% 
                                                                   names(data)))
        col[which(sampleID %in% names(data))] <- "firebrick"
    }

    opar <- par(cex=1.6)
    on.exit(par(opar))

    plot(seq_along(data), data, col = col, pch=20,
        xlab="sample rank", ylab=dist, las=1,
        main=paste("Distribution of:", dist)
    )
}

# psiType="psi3"
# dist     <- "psi3"
# sampleID <- "sample1"
# idx      <- 48


#'
#' Present the FraseR results as shiny app
#'
#' @return shiny application
#' @examples
#'     fds <- createTestFraseRDataSet()
#'
#'     # for interactive sessions
#'     myShinyApp <- FraseRShinyApp(fds)
#'
#'     # for running a shiny application as server
#'     myShinyAppObj <- FraseRShinyApp(fds, server=TRUE)
#'
#' @noRd
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
