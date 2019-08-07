#' @export

plotSeq3Dims <- function(tforms){
  invisible(lapply(tforms, function(x){
    if (!is.null(x$overlay)){
      display(x$overlay)
    }
  }))
}
