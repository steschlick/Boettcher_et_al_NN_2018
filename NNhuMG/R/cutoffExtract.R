#' @title Extract Gate Boundaries from GatingSet
#' @description convenience function 
#' @param x matrix with gate vertices
#' @param yChnl channel against which target marker has been plotted in Cytobank 
#' @param xyChnl optional channel against which yChnl was plotted in Cytobank
#' @return cutoff value
#' @export 
cutoffExtract <- function(x, yChnl, xyChnl) {
  cols <- colnames(x)
  if (all(cols %in% c(yChnl, xyChnl))) {
    co <- min(x[, cols == yChnl])
    names(co) <- yChnl
  } else {
    co <- min(x[, cols != yChnl])
    names(co) <- cols[cols != yChnl]
  }
  co
}
