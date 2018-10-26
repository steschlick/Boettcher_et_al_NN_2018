#' @title Scale Weights for Plotting
#' @description convenience function
#' @param w numeric vector
#' @param th numeric threshold
#' @return numeric vector
#' @details Sets values below threshold to zero. Scales weights to sum up to the number of observations. 
#' @export 
threshold2zero <- function(w, th) {
  w <- w/th
  w <- w - 1
  w[w < 0] <- 0
  w <- w * length(w)/sum(w)
  w
}
