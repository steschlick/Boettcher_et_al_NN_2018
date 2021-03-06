#' @title Robust Effect Size
#' @description computes a Hodges-Lehmann estimator-based effect size
#' @param x numeric vector of length m
#' @param y numeric vector of length n
#' @param th numeric threshold value for the scale of the two distributions 
#' to avoid unintentional large \eqn{\Delta} values when expression levels are near zero
#' @return \eqn{\Delta} numeric value
#' @details \eqn{\Delta} is defined as the Hodges-Lehmann two-sample estimator scaled by the median absolute deviation about the 
#'          two-sampleHodges-Lehmann estimator. Currently, only input vectors with maximum length 1024 are allowed for feasible computation.
#'          If their size exceeds this limit, data can be downsampled using \code{\link{exprs2con}}.
#' @examples 
#' \dontrun{
#' # without outliers, compare to Cohen's d
#' Cohen <- function(x,y) {
#'   n <- length(x)-1
#'   m <- length(y)-1
#'   d  <- mean(x)-mean(y)
#'   pv <- (n*var(x)+m*var(y))/(n+m)
#'   d/sqrt(pv)           
#' }
#' x <- rnorm(512, 10, 2)             
#' y <- rnorm(512, 14, 2)
#' Cohen(x, y)
#' robDelta(x, y)
#' 
#' # with outliers
#' x <- c(rnorm(480, 10, 2),  rnorm(32, 21, 2))             
#' Cohen(x, y)
#' robDelta(x, y)
#' }
#' @references 
#'   P. Rousseeuw and C. Croux (1992), Explicit Scale Estimators with High Breakdown point, L1-Statistical Analysis and Related Methods
#' @seealso 
#'  \code{\link[stats]{median}}
#'  \code{\link[stats]{mad}}
#' @rdname robDelta
#' @export 
#' @importFrom stats median mad
robDelta <- function(x, y, th = 0) {
  if (length(x) > 1024L || length(y) > 1024L) {
   stop("too many observations to compute delta - use expr2con() for downsampling") 	
  }
  out <- outer(x, y, "-")
  m <- stats::median(out)  # HL
  return(sqrt(2)*m/max(stats::mad(out), th))
}
