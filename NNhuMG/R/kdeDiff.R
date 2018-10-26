#' @title Positive Density Difference
#' @description subtracts two densities and returns the positive difference
#' @param x matrix of data values
#' @param fhat object returned by \code{\link{kde}}
#' @param binned flag for binned estimation, Default: TRUE
#' @param bgridsize vector of binning grid sizes
#' @param H bandwidth matrix
#' @param w0 vector of weights for the minuend density, Default is a vector of all ones.
#' @param w1 vector of weights for the subtrahend density
#' @param compute.cont flag for computing 0.1 to 0.99 probability contour levels, Default: FALSE
#' @return object of class \code{\link{kde}} with the positive difference of two density estimates in field \code{fhat}
#' @seealso 
#'  \code{\link[ks]{kde}}
#' @rdname kdeDiff
#' @export 
#' @importFrom ks kde
kdeDiff <- function(x, fhat, binned = TRUE, bgridsize, H, w0, w1, compute.cont = FALSE) {
  
  if (!missing(fhat)) {
    f0 <- fhat
  } else {
    if (missing(w0)) {
      f0 <- ks::kde(x = x, binned = binned, bgridsize = bgridsize, H = H)
    } else f0 <- ks::kde(x = x, binned = binned, bgridsize = bgridsize, H = H, w = w0)
  }
  
  if (!missing(w1)) {
    if (missing(x)) {
      x <- fhat$x
      binned <- fhat$binned
      bgridsize <- sapply(fhat$eval.points, length)
      H <- fhat$H
    }
    f1 <- ks::kde(x = x, binned = binned, bgridsize = bgridsize, H = H, w = w1)
    f0$estimate <- f1$estimate - f0$estimate
    f0$estimate[f0$estimate < 0] <- 0
    f0
  } else stop("w1 is required")
  
}
