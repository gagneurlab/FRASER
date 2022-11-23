#'
#' Imported and reexported generic functions
#' Mainly taken from OUTRIDER and generics
#' 
#' @noRd
NULL


#' @importFrom OUTRIDER aberrant
#' @noRd
#' @export
OUTRIDER::aberrant

#' @importFrom OUTRIDER filterExpression
#' @noRd
#' @export
OUTRIDER::filterExpression

#' @importFrom generics fit
#' @noRd
#' @export
generics::fit

#' @importFrom OUTRIDER filterExpression
#' @noRd
#' @export
OUTRIDER::filterExpression

#' @importFrom OUTRIDER plotAberrantPerSample
#' @noRd
#' @export
OUTRIDER::plotAberrantPerSample

#' @importFrom OUTRIDER plotCountCorHeatmap
#' @noRd
#' @export
OUTRIDER::plotCountCorHeatmap

#' @importFrom OUTRIDER plotEncDimSearch
#' @noRd
#' @export
OUTRIDER::plotEncDimSearch

#' @importFrom OUTRIDER plotQQ
#' @noRd
#' @export
OUTRIDER::plotQQ

#' @importFrom OUTRIDER plotVolcano
#' @noRd
#' @export
OUTRIDER::plotVolcano

#' @importFrom OUTRIDER results
#' @noRd
#' @export
OUTRIDER::results



#'
#' FRASER specific functions defined as generics
#' 
#' @noRd
NULL


#' @rdname fds-methods
#' @export
setGeneric("samples", function(object) standardGeneric("samples"))

#' @rdname fds-methods
#' @export
setGeneric("samples<-", signature = "object", 
        function(object, value) standardGeneric("samples<-"))

#' @rdname fds-methods
#' @export
setGeneric("condition", function(object) standardGeneric("condition"))

#' @rdname fds-methods
#' @export
setGeneric("condition<-", signature = "object", 
        function(object, value) standardGeneric("condition<-"))

#' @rdname fds-methods
#' @export
setGeneric("bamFile", function(object) standardGeneric("bamFile"))

#' @rdname fds-methods
#' @export
setGeneric("bamFile<-", signature = "object",
        function(object, value) standardGeneric("bamFile<-"))

#' @rdname fds-methods
#' @export
setGeneric("name", function(object) standardGeneric("name"))

#' @rdname fds-methods
#' @export
setGeneric("name<-", signature = "object", 
        function(object, value) standardGeneric("name<-"))

#' @rdname fds-methods
#' @export
setGeneric("strandSpecific", function(object) standardGeneric("strandSpecific"))

#' @rdname fds-methods
#' @export
setGeneric("strandSpecific<-",  signature = "object", 
        function(object, value) standardGeneric("strandSpecific<-"))

#' @rdname fds-methods
#' @export
setGeneric("pairedEnd", function(object) standardGeneric("pairedEnd"))

#' @rdname fds-methods
#' @export
setGeneric("pairedEnd<-", signature = "object", 
        function(object, value) standardGeneric("pairedEnd<-"))

#' @rdname fds-methods
#' @export
setGeneric("workingDir", function(object) standardGeneric("workingDir"))

#' @rdname fds-methods
#' @export
setGeneric("workingDir<-", signature = "object", 
        function(object, value) standardGeneric("workingDir<-"))

#' @rdname fds-methods
#' @export
setGeneric("scanBamParam", function(object) standardGeneric("scanBamParam"))

#' @rdname fds-methods
#' @export
setGeneric("scanBamParam<-", signature = "object", 
        function(object, value) standardGeneric("scanBamParam<-"))

#' @rdname fds-methods
#' @export
setGeneric("nonSplicedReads", 
        function(object) standardGeneric("nonSplicedReads"))

#' @rdname fds-methods
#' @export
setGeneric("nonSplicedReads<-", signature = "object", 
        function(object, value) standardGeneric("nonSplicedReads<-"))

#' @rdname plotFunctions
#' @export
setGeneric("plotManhattan", function(object, ...) 
    standardGeneric("plotManhattan"))

