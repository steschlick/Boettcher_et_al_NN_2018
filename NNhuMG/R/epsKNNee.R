#' @title Estimate EPS
#' @description robust knee detection to estimate a suitable value for dbscan eps neighborhood
#' @param x matrix of observations
#' @param k number of nearest neighbors
#' @param plotit flag to plot the knee, Default: FALSE
#' @return eps, size of the epsilon neighborhood
#' @examples 
#' \dontrun{
#' # example data from fpc
#' set.seed(665544)
#' n <- 100
#' x <- cbind(
#'   x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
#'   y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
#'   )
#' epsKNNee(x, k = 3, plotit=TRUE)
#' }
#' @seealso 
#'  \code{\link[dbscan]{kNNdist}}
#'  \code{\link[dbscan]{kNNdistplot}}
#'  \code{\link[robust]{lmRob}}
#' @rdname epsKNNee
#' @export 
#' @importFrom dbscan kNNdist
#' @importFrom robust lmRob
#' @importFrom graphics plot abline
epsKNNee <- function(x, k, plotit = FALSE) {
  h <- dbscan::kNNdist(x, k = k)
  h <- sort(h)
  n <- length(h)
  h.i <- seq(1, n)
  p <- which.max(sapply(h.i, function(x) {
    sqrt((h[x] - h[1])^2 + (x - 1)^2) + sqrt((h[n] - h[x])^2 + (n - x)^2)
  }))
  p <- min(p, n - 5)
  p <- max(p, 7)
  df1 <- data.frame(x = h.i[1:(p - 1)], y = h[1:(p - 1)])
  df2 <- data.frame(x = h.i[p:n], y = h[p:n])
  coef1 <- robust::lmRob(y ~ x, data = df1)$coefficients
  coef2 <- robust::lmRob(y ~ x, data = df2)$coefficients
  eps <- unname((coef1[1]/coef1[2] - coef2[1]/coef2[2])/(1/coef1[2] - 1/coef2[2]))
  if (plotit) {
    graphics::plot(h, pch = ".", ylab = "distance")
    graphics::abline(coef1, col = 2)
    graphics::abline(coef2, col = 2)
    graphics::abline(h = eps, col = 3)
  }
  eps
}
