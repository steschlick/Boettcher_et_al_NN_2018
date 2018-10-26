#' @title Scale Weights for Density Difference
#' @description a ramp function
#' @param w numeric vector, either effect size or -log transformed \emph{p}-values
#' @param th numeric threshold, e.g. -log transformed level of significance
#' @return numeric vector of weights summing up to the number of observations
#' @details The function piece-wise linear transforms input values such that they are 1 if equal to \code{th} and 
#'       their sum equals the number of observations. 
#' @examples 
#' \dontrun{
#'  # vector of p-values
#'  p <- runif(512)
#'  w <- threshold2one(-log2(p), -log2(0.05))
#'  sum(w)
#'  plot(sort(threshold2one(1:100, 98)))
#'  plot(sort(threshold2one(1:100, 50)))
#'  plot(sort(threshold2one(1:100, 2)))
#' }
#' @rdname threshold2one
#' @export 
threshold2one <- function(w, th) {
  n <- length(w)
  w <- w/th
  w.neg <- w < 1
  w.pos <- w >= 1
  n.neg <- sum(w.neg)
  n.pos <- sum(w.pos)
  swit <- n.pos <= n/2
  if (swit) {
    w <- w - 1
    w[w.pos] <- n.pos * w[w.pos]/(2 * sum(w[w.pos]))
    w[w.neg] <- w[w.neg] * sum(w[w.pos])/sum(abs(w[w.neg]))
    w <- w + 1
  } else {
    w[w.pos] <- w[w.pos] - 1
    target.pos <- n - sum(abs(w[w.neg])) - n.pos
    w[w.pos] <- target.pos * w[w.pos]/sum(abs(w[w.pos]))
    w[w.pos] <- w[w.pos] + 1
  }
  w
}