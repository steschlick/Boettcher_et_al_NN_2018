#' @title Peak detection
#' @description Finds peaks in a 2D kernel density estimate
#' @param fhat object returned by \code{\link{kde}} or \code{\link{kdeDiff}}
#' @param tol threshold to avoid detection of noise, Default: 5e-08
#' @return peak points matrix with local densities in 3rd column  
#' @details This function identifies local maxima in a kernel density estimate. 
#' Code has been adaptated from the ACCENSE R-implementation accessible under \url{http://www.cellaccense.com/oldver.html}. 
#' @examples 
#' \dontrun{
#' library(ks)
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
#' fhat <- kde(x)
#' pts <- PeakFind(fhat)
#' plot(fhat)
#' points(pts[ , 1:2], pch="*")
#' }
#' @references 
#'   K. Shekhar, P. Brodin, M.M. Davis, A.K. Chakraborty (2014), Automatic classification of 
#'   cellular expression by nonlinear stochastic embedding (ACCENSE), 
#'   Proceedings of the National Academy of Sciences, \url{http://www.cellaccense.com/oldver.html}, 
#' 
#'   E.R. Davies (1997), Machine Vision, volume 609, Academic press New York.
#' @rdname PeakFind
#' @export
#' @keywords internal 
PeakFind <- function(fhat, tol = 5e-08) {
# code adapted from http://cytof.scilifelab.se/accenseweb/download/release_abstract.wenjian.861cb69fbbd18a0a.525f696d706c656d656e746174696f6e2e7a6970.zip

  d <- fhat$estimate
  krx <- kronecker(matrix(1, 1, dim(d)[2]), matrix(fhat$eval.points[[1]], ncol = 1))
  kry <- kronecker(matrix(1, dim(d)[1], 1), matrix(fhat$eval.points[[2]], nrow = 1))  
  edg <- 1  # needs to be min 1
  threshold <- max(c(min(apply(d, 1, max)), min(apply(t(d), 1, max))), tol)
  # Apply threshold
  indicator <- (d > threshold) * 1  #Indicator matrix
  d <- d * indicator
  if (sum(d) != 0) {
     rows_cols <- which(d[(1 + edg):(dim(d)[1] - edg), (1 + edg):(dim(d)[2] - 
      edg)] > 0, arr.ind = TRUE)
    x <- rows_cols[, 1]
    y <- rows_cols[, 2]
    cent <- c(0, 0, 0)  # Initialize output
    x <- x + edg
    y <- y + edg
    for (j in 1:length(y)) {
      if (d[x[j], y[j]] >= d[x[j] - 1, y[j] - 1] && d[x[j], y[j]] > d[x[j] - 
        1, y[j]] && d[x[j], y[j]] >= d[x[j] - 1, y[j] + 1] && d[x[j], y[j]] > 
        d[x[j], y[j] - 1] && d[x[j], y[j]] > d[x[j], y[j] + 1] && d[x[j], 
        y[j]] >= d[x[j] + 1, y[j] - 1] && d[x[j], y[j]] > d[x[j] + 1, y[j]] && 
        d[x[j], y[j]] >= d[x[j] + 1, y[j] + 1]) {
        cent <- rbind(cent, c(krx[x[j], y[j]], kry[x[j], y[j]], fhat$estimate[x[j], 
          y[j]]))
      }
    }
    # Remove first row which is a dummy
    cent <- matrix(cent[-1, ], ncol = 3)
    colnames(cent) <- c("x", "y", "d")
    rownames(cent) <- 1:nrow(cent)
    
    return(cent)
  } else return(NA)
  
}
