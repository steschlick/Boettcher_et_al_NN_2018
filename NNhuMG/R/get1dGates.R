#' @title 1D Gates
#' @description creates 1d flowCore filters for 1, 2, or 3 cutoffs and a given marker
#' @param chnl channel name
#' @param mrkr marker name
#' @param cutoffs numeric vector of cutoff-values
#' @param boundaries a lower and an upper limit of the expression data, e.g. the range
#' @return A list containing the flowCore filter objects
#' @details The function creates a set of all possible 1d-\code{\link{rectangleGate}s} for a given number of cutoffs and assigns filterId's according 
#'          to the phenotyping syntax specified in \code{\link{phenoCutoff}}. Filters of different markers may be combined to create multidimensional 
#'          'boolean' gates to target specific phenotypes.
#' @examples 
#' \dontrun{
#'  #EXAMPLE1
#'  # in flowCore
#'  cd3 <- get1dGates("Eu151Di", "CD3", 4, c(0, 8))
#'  cd4 <- get1dGates("Yb173Di", "CD4", c(2, 5), c(0, 8))
#'  cd3$`CD3-` * cd4$`CD4+-`
#'  # add to GatingSet
#'  exp <- "/CD3- & /CD4+-"
#'  bf <- boolFilter(exp, filterId = "Mono")
#'  # add(gs, cd3$`CD3-`)
#'  # add(gs, cd4$`CD4+-`)
#'  # add(gs, bf)
#' }
#' @seealso 
#'  \code{\link[flowCore]{rectangleGate-class}}
#'  \code{\link{boolFilter}}
#' @rdname get1dGates
#' @export 
#' @importFrom flowCore rectangleGate
get1dGates <- function(chnl, mrkr, cutoffs, boundaries) {
  n <- length(cutoffs)  # cuts should be sorted already
  if (n < 1 || n > 3) {
    NULL
  } else {
    switch(n, 
    stats::setNames(
    list(flowCore::rectangleGate(filterId = paste0(mrkr, "-"), .gate = matrix(c(boundaries[1], cutoffs[1]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "+"), .gate = matrix(c(cutoffs[1], boundaries[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))) 
    ), paste(mrkr, c("-", "+"), sep="")), 
    stats::setNames(
    list(flowCore::rectangleGate(filterId = paste0(mrkr, "--"), .gate = matrix(c(boundaries[1], cutoffs[1]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "-"), .gate = matrix(c(boundaries[1], cutoffs[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "+-"), .gate = matrix(c(cutoffs[1], cutoffs[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "+"), .gate = matrix(c(cutoffs[1], boundaries[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "++"), .gate = matrix(c(cutoffs[2], boundaries[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))) 
    ), paste(mrkr, c("--", "-", "+-", "+", "++"), sep="")),  
    stats::setNames(
    list(flowCore::rectangleGate(filterId = paste0(mrkr, "---"), .gate = matrix(c(boundaries[1], cutoffs[1]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "--"), .gate = matrix(c(boundaries[1], cutoffs[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "-"), .gate = matrix(c(boundaries[1], cutoffs[3]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "+--"), .gate = matrix(c(cutoffs[1], cutoffs[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "+-"), .gate = matrix(c(cutoffs[1], cutoffs[3]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "++-"), .gate = matrix(c(cutoffs[2], cutoffs[3]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "+"), .gate = matrix(c(cutoffs[1], boundaries[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "++"), .gate = matrix(c(cutoffs[2], boundaries[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
         flowCore::rectangleGate(filterId = paste0(mrkr, "+++"), .gate = matrix(c(cutoffs[3], boundaries[2]), ncol = 1, dimnames = list(c("min", "max"), c(chnl)))), 
    ), paste(mrkr, c("---", "--", "-", "+--", "+-", "++-", "+", "++", "+++"), sep="")))
  }
}






