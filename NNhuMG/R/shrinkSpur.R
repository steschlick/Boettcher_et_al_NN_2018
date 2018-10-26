#' @title Shrink Spurious Contours
#' @description removes outer contours that have been reported due to presence of spurious peaks
#' @param con object returned by \code{\link{contourLoc}}
#' @param selected indices of contours corresponding to rows in \code{con$res}, useful to manually select or deselect contours
#' @param mul if TRUE, \code{minBins} is multiplied by the number of inner spurious peaks (as there might be multiple) to more 
#'        progressively shrink towards single inner contours, Default: FALSE
#' @param minBins minimum tolerated difference in the numbers of observations of an outer and its single inner contour.
#' @param peaks2merge may also be lowered in order to remove outmost contours
#' @return  \code{con} list object with added field:
#' \describe{
#'   \item{\code{shrunk}}{indices of kept contours}
#' }
#' @details \code{\link{contourLoc}} might have returned contours that surround a single inner contour (if \code{peaks2merge > 1}), 
#'        due to the detection of spurious peaks. These contours need to be removed in order to have a sufficient number of observations
#'        (specified by \code{minBins}) for comparing marker expression levels of contour-gated subsets as performed by \code{\link{contourMerge}}.
#'        Contours are labeled by tupels of peak indices, so they can be visually identified in a plot (see \code{\link{plotContours}}). Note that 
#'        currently no print, show, or summary functions are implemented.
#' @examples 
#' \dontrun{
#'  set.seed(123) 
#'  x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2, 
#'                    dimnames=list(NULL, c("A", "B"))),
#'             matrix(rnorm(36, mean = 1, sd = 0.3), ncol = 2),
#'             matrix(rnorm(100, mean = -1, sd = 0.3), ncol = 2))
#'  BW <- chooseBandWidth(x) 
#'  
#'  # get density estimate
#'  fhat <- BW$fhats[[ BW$Hindex[1] ]]
#'  
#'  # get contour list
#'  cs <- contourLoc(fhat=fhat, offset = 0.01, n.levels = 20, peaks2merge = 5, minBins = 6,
#'                   keep.level.by = "max.bin.num", plotit = FALSE)
#'  # spurious peak no. 5 
#'  plotContours(cs, cex.labels=1.5, fhat.levels=T, labels="peak.ids")
#'  
#'  # remove outer contours containing too few observations
#'  css <- shrinkSpur(cs) 
#'  plotContours(css, fhat.levels=T, labels="peak.ids")
#' }
#' @rdname shrinkSpur
#' @export 
shrinkSpur <- function(con, selected, mul = FALSE, minBins, peaks2merge) {
  
  if (missing(minBins))
    minBins <- con$minBins else con$minBins <- minBins
  
  if (missing(selected)) 
    {
      selected <- con$update
    }  # else {
  # if (is.numeric(selected)) selected <- seq_len(nrow(con$res)) %in% selected }
  
  obj <- con$res[selected, ]
  
  # adjacency matrix
  A <- contourAdjacency(obj$peaks)
  
  # contours to keep
  keep <- logical(ncol(A))
  
  # contours w/o inner contours
  no.in <- apply(A, 2, function(x) !any(x != 0))
  # keep innermost contours by default
  keep[no.in] <- TRUE
  
  # remove outer contours having single inner contour and difference in bin number
  # < minBins, since these cannot be evaluated with regard to their phenotypic
  # differences
  topo <- sort(unique(A[1:length(A)]))
  if (!missing(peaks2merge)) 
    topo <- topo[topo <= peaks2merge]
  topo <- rev(topo[-1])
  
  for (t in topo) {
    co <- which(apply(A, 2, function(x) any(x == t)))
    for (c in co) {
      inner <- which(as.logical(A[, c]))
      inner.n <- length(inner)
      if ((inner.n == 1) && ((obj$num.bins[[c]] - obj$num.bins[[inner]]) < 
        (ifelse(mul, t - length(obj$peaks[[inner]]), 1) * minBins))) 
        next else keep[c] <- TRUE
    }
  }
  
  # update object
  con$update <- selected[keep]
  con$shrunk <- selected[keep]
  
  return(con)
}
