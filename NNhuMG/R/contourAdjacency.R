## @title helper function
## @description computes adjacency matrix for contours
## @param peaks a list whose elements are vectors of peak indices for each contour
## @return adjacency matrix of where 1s indicate an 'inner adjacency'
## @details We use peak indices to represent contours as well as to define their membership with regard to inner and outer. This is possible 
##          since \code{\link{contourLoc}} returns contours that never share a common tuple of peaks. 
#' @keywords internal
contourAdjacency <- function(peaks) {
  A <- diag(0, length(peaks))
  dimnames(A) <- list(names(peaks), names(peaks))
  for (i in 1:length(peaks)) {
    l <- length(peaks[[i]])
    if (l == 1) 
      next
    g <- sapply(peaks, function(x) all(x %in% peaks[[i]]))
    s <- sum(g)
    if (s > 1) {
      m <- sapply(1:s, function(j) {
        sum(sapply(peaks[g][-j], function(y) {
          all(peaks[g][[j]] %in% y)
        })) == 1
      })
      g[g][!m] <- FALSE
    }
    g[i] <- FALSE
    A[g, i] <- l
  }
  A
}
