#' @title \emph{comb}inatorial \emph{Comp}arisons
#' @description compares marker expression levels between each pair of selected subsets
#' @param res data.frame with column \code{bin.idx}, returned by \code{\link{contourLoc}}
#' @param sel row indices to select in res, i.e. the subsets 
#' @param data matrix containing marker expression values
#' @param delta numeric threshold for robust effect size \eqn{\Delta}, Default: NA
#' @param beta numeric threshold for non-overlapping quantiles, Default: 0.25
#' @param threshold numeric for distance between two distributions, applied at the same scale 
#' of marker expression, this also avoids unintentional large \eqn{\Delta} values when expression levels are near zero, 
#' passed to \code{\link{robDelta}} and \code{\link{phenoQuant}}, Default: 0.5
#' @return A list containing elements:
#' \describe{
#'   \item{idx}{list of contour/subset indices for pairwise comparisons}
#'   \item{keep}{logical matrix indicating the results with markers in columns for which subsets have been compared}
#'   \item{delta}{numeric matrix with computed signed \eqn{\Delta} values, same dimensions as \code{keep}}
#'   \item{beta}{numeric matrix with signs indicating non-overlap and direction, same dimensions as \code{keep}}
#' }
#' @details This helper function compares each pair of selected subsets (or contours) for phenotypic similarity for  
#'  each marker found in \code{data} by means of two dicision boundaries: quantile overlap and/or a minimally   
#'  required robust effect size \eqn{\Delta}. Since the former does not provide information on how much two subsets 
#'  differ with regard to a given marker, the latter may be used to estimate the importance of each marker to 
#'  differentiate between two subsets (as in \code{\link{phenoReduce}}). As possible results, markers that do not
#'  provide distinction between any two subsets are being excluded from downstream analysis or two phenotypically 
#'  similar subsets might be merged (as in \code{\link{contourMerge}}).
#' @seealso 
#'  \code{\link{contourLoc}}
#'  \code{\link{contourMerge}}
#'  \code{\link{robDelta}}
#'  \code{\link{exprs2con}}
#'  \code{\link{findCutoffs}}
#'  \code{\link{phenoQuant}}
#' @importFrom utils combn
#' @importFrom stats setNames
#' @keywords internal
combCompare <- function(res, sel, data, delta = NA, beta = 0.25, threshold = 0.5) {
  # gsub('_|-|\\+|&|/|\\|| \\(v\\)', '', colnames(data))
  
  use.delta <- !is.na(delta)
  use.beta <- !is.na(beta)
  if (!use.delta && !use.beta) 
    stop("delta and/or beta need to be provided")
  
  comb.idx <- sapply(sel, function(x) {
    utils::combn(x, 2, simplify = FALSE)
  }, simplify = FALSE)
  
  comb.idx <- unlist(comb.idx, recursive = F, use.names = F)
  sel.nms <- rownames(res)
  names(comb.idx) <- sapply(comb.idx, function(x) {
    paste(sel.nms[x[1]], sel.nms[x[2]], sep = " vs ")
  })
  
  if (use.delta) {
    delta.mat <- t(sapply(comb.idx, function(i) {
      sapply(1:ncol(data), function(j) {
          robDelta(x = data[res$bin.idx[[i[1]]], j], 
                y = data[res$bin.idx[[i[2]]], j], th = threshold)
      })
    }, simplify = TRUE))
    colnames(delta.mat) <- colnames(data)
  } else delta.mat <- NULL
  
  if (use.beta) {
    beta.mat <- t(sapply(comb.idx, function(i) {
      sapply(1:ncol(data), function(j) {
        phenoQuant(x = data[res$bin.idx[[i[1]]], j], y = data[res$bin.idx[[i[2]]], 
          j], beta = beta, th = threshold)
      })
    }, simplify = T))
    colnames(beta.mat) <- colnames(data)
  } else beta.mat <- NULL
  
  if (use.delta && use.beta) {
    keep.mat <- (abs(delta.mat) >= delta) & (abs(beta.mat) == 1)
  } else {
    if (use.delta) {
      keep.mat <- abs(delta.mat) >= delta
    } else keep.mat <- abs(beta.mat) == 1
  }
  # 
  return(stats::setNames(list(comb.idx, keep.mat, delta.mat, beta.mat), c("idx", "keep", 
    "delta", "beta")))
  
}
