#' @title Merge Similar Subsets
#' @description evaluates phenotypic similarity of contour-gated subsets. An outer contour is removed if it contains 
#'        inner subsets that have different phenotype(s). Otherwise, inner subsets are "merged" with the outer.  
#' @param con object returned by \code{\link{shrinkSpur}}
#' @param selected indices of contours to be compared, useful to manually deselect contours
#' @param keep.inner if TRUE, inner contours are not merged even though there are similar (to the outer), Default: FALSE
#' @param keep.outer if TRUE, outer contours are not removed even though there are different from inner, Default: FALSE
#' @param minBins minimum difference in numbers of observations wihin an outer and inner contours to consider the exclusion (setdifference)
#'        of inner and outer for comparisons (see \code{\link{shrinkSpur}}).
#' @param data matrix containing marker expression values
#' @param delta numeric threshold for robust effect size \eqn{\Delta}, Default: NA
#' @param beta numeric threshold for non-overlapping quantiles, Default: 0.25
#' @param threshold numeric for distance (at asinh scale) between two distributions to avoid unintentional 
#'        large \eqn{\Delta} values when expression levels are near zero, Default: 0.5
#' @return \code{con} list object with added fields:
#' \describe{
#'   \item{\code{merged}}{indices of kept contours}
#'   \item{\code{data}}{expression matrix}
#'   \item{\code{comb}}{\describe{
#'                        \item{\code{sel}}{data.frame containing indices for each set of inner and outer contours}
#'                        \item{\code{res}}{list of results for each set of inner and outer comparisons as returned by \code{\link{combCompare}}}
#'                      }}
#' }
#' @details Contour-gated subsets are compared for phenotypic similarity with regard to marker expression provided in \code{data}. 
#'        They will be considered dissimilar if, for any marker, their symmetric quantiles (as specified by \code{beta}) overlap and/or  
#'        the expression exceeds the limit for robust effect size \eqn{\Delta}. Comparisons are performed using \code{\link{combCompare}} for  
#'        all outmost contours and their exclusion (i.e. the "root") as well as for each set of outer and inner contours. The results of the
#'        final comparisons between dissimilar subsets are returned with the contour list object and may be inspected via the 
#'        high-level plotting function \code{\link{plotPheno}}. Note that currently no print, show, or summary functions are implemented.
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
#'                   keep.level.by = "max.bin.num", plotit = TRUE)
#'  
#'  # remove contours containing too few observations
#'  css <- shrinkSpur(cs) 
#'  
#'  # compare and keep outer subsets
#'  csm <- contourMerge(css, keep.outer=T, keep.inner=F, data=x, 
#'                      delta=2, beta=0.15, threshold=0.2)
#'  plotContours(csm, col="grey", pch="+", cex=0.5, lwd=2.5, fhat.levels=T, labels="peak.ids")
#'  
#'  # available subsets
#'  rownames(csm$res)
#'  
#'  # sets of outer and inner contours
#'  csm$comb$sel
#'  
#'  # compare within outer 5:
#'  plotContours(csm, comp.sel=c(1, 4, 8), col="grey", pch="+", cex=0.5, lwd=2.5, fhat.levels=T, 
#'               labels="peak.ids", filwd=.5, density=50)
#'  plotPheno(csm, comp.sel=c(1, 4, 8)) # inner differ phenotypically
#'  
#'  # outer removed
#'  csm <- contourMerge(css, keep.outer=F, keep.inner=F, data=x, 
#'                      delta=2, beta=0.15, threshold=0.2)
#'  plotContours(csm, col="grey", pch="+", cex=0.5, lwd=2.5, fhat.levels=T, labels="peak.ids")
#'  
#'  # show effect size as well
#'  plotPheno(csm, comp.sel=T, ask=F) 
#' }
#' @seealso 
#'  \code{\link{combCompare}}
#'  \code{\link{plotPheno}}
#' @rdname contourMerge
#' @export 
#' @importFrom stats setNames
contourMerge <- function(con, selected, keep.inner = FALSE, keep.outer = FALSE, minBins, 
  data, delta = NA, beta = 0.25, threshold = 0.5) {
  
  # some checks and arg updates
  if (is.null(con$shrunk)) 
    stop("use shrinkSpur() first to remove contours enclosing spurious peaks")
  
  use.delta <- !is.na(delta)
  use.beta <- !is.na(beta)
  
  if (!use.delta && !use.beta) 
    stop("a value for delta and/or beta needs to be provided")
  con$delta <- delta
  con$beta <- beta
  con$threshold <- threshold
  
  if (missing(selected)) {
    selected <- con$shrunk  # shrunken
  } else {
    # here test for correct indices TODO: provide only if used in interactive mode?
    # implement update function for interactive selection!
    if (!all(selected %in% con$shrunk)) 
      stop("wrong indices selected")
  }
  
  if (missing(minBins))
    minBins <- con$minBins else con$minBins <- minBins

  
  if (!(keep.inner && keep.outer)) {
    
    # first, remove inner spurious contours, i.e. that won't be distinguishable
    # phenotypically using combCompare() based on a set of markers. In case they are
    # different, remove the outer contour.
    
    # temporary contour res list
    obj <- con$res[selected, ]
    
    A <- contourAdjacency(obj$peaks)
    
    keep <- logical(ncol(A))
    
    # innermost contours
    no.in <- apply(A, 2, function(x) !any(x != 0))
    
    # flag innermost contours
    keep[no.in] <- TRUE
    
    # contains inner contour(s)
    co.in <- apply(A, 2, function(x) any(x != 0))
    
    # from inner to outer
    topo <- sort(unique(A[1:length(A)]))
    topo <- topo[-1]
    
    for (t in topo) {
      co <- which(apply(A, 2, function(x) any(x == t)))
      for (c in co) {
        inner <- which(as.logical(A[, c]))
        if (any(!keep[inner])) {
          keep[c] <- FALSE
          # next
        } else {
          
          co.bins <- obj$bin.idx[c]
          inner.bins <- obj$bin.idx[inner]
          incl.bins <- unique(unname(unlist(inner.bins)))
          excl.bins <- list(co.bins[[1]][!co.bins[[1]] %in% incl.bins])
          
          tmp.res <- list()
          
          # don't actually need outer contour, i.e. obj$bin.idx[c], but keep for plotting
          # TODO: if length(inner) == 1 and inner/outerexcl differ, then keep the excl
          if ((length(inner) == 1) | (length(excl.bins[[1]]) >= minBins)) {
          excl.peaks <- list(obj$peaks[[c]])
          # excl.peaks <- list(setdiff(obj$peaks[[c]],
          # unique(unname(unlist(obj$peaks[inner])))))
          names(excl.bins) <- paste(c("*", excl.peaks[[1]]), collapse = "\u00B7")
          tmp.res$bin.idx <- c(obj$bin.idx[inner], excl.bins)
          
          } else {
          tmp.res$bin.idx <- c(obj$bin.idx[inner])
          }
          
          keep.mat <- combCompare(res = as.data.frame(do.call("cbind", tmp.res)), 
          data = data, sel = list(tmp = seq_len(length(tmp.res$bin.idx))), 
          delta = delta, beta = beta, threshold = threshold)$keep
          marker.keep <- apply(keep.mat, 2, any)
          
          if (any(marker.keep)) {
          keep[c] <- keep.outer
          keep[inner] <- TRUE
          } else {
          keep[c] <- TRUE
          keep[inner] <- keep.inner
          }
        }
      }
    }
    
    selected <- selected[keep]
    
  } else message("object will be returned w/o effect when both keep.inner and keep.outer are TRUE")
  
  # second, update the con res object for comparisons of merged contours and
  # visualization update object
  con$update <- selected
  con$merged <- selected
  
  # adjacency for final comparisons
  A <- contourAdjacency(con$res[selected, ]$peaks)
  
  # holds res indices of contours to be compared
  sel.comb <- list()
  
  # setup 'root' comparison of all outmost contours
  no.out <- apply(A, 1, function(x) !any(x != 0))
  root <- selected[no.out]
  # collect the exclusion of all outmost contours as additional reference against
  # which separating markers/cutoffs need to be found as well if (nrow(data) !=
  # nrow(con$fhat$x)) stop('row numbers of parameter and tsne data do not match')
  all.bins <- seq_len(nrow(data))
  inner.bins <- con$res$bin.idx[root]
  incl.bins <- unique(unname(unlist(inner.bins)))
  excl.bins <- all.bins[!all.bins %in% incl.bins]
  excl.peaks <- setdiff(seq_len(nrow(con$peaks)), unique(unname(unlist(con$res$peaks[root]))))
  if (!length(excl.peaks)) 
    excl.peaks <- 0
  excl.names <- paste(c("*\u00B7root", excl.peaks), collapse = "\u00B7")
  
  # get index of exclusion if already present
  excl.idx <- which(rownames(con$res) %in% excl.names)
  
  if (!length(excl.idx)) {
    # add to res list
    con$res <- rbind(con$res, do.call("rbind", stats::setNames(list(tmp = list(peaks = excl.peaks, 
      p.value = NA_real_, boundary = list(x = NA_real_, y = NA_real_), bin.idx = excl.bins, 
      num.bins = length(excl.bins), eff.bins = NA_real_)), excl.names)))
    excl.idx <- nrow(con$res)
  }
  # and keep indices
  sel.comb[["root"]] <- list(outer = NA_real_, inner = c(root, excl.idx))
  
  # same for all contours that contain inner contour(s)
  co.in <- apply(A, 2, function(x) any(x != 0))
  
  co <- which(co.in)
  
  for (c in co) {
    
    inner <- selected[which(as.logical(A[, c]))]
    outer <- selected[c]
    
    # use minBins here!
    co.bins <- con$res$bin.idx[outer]
    inner.bins <- con$res$bin.idx[inner]
    incl.bins <- unique(unname(unlist(inner.bins)))
    excl.bins <- co.bins[[1]][!co.bins[[1]] %in% incl.bins]
    
    if ((length(inner) == 1) | (length(excl.bins) >= minBins)) {
      
      # collect the exclusion within outer contour as reference excl.peaks <-
      # setdiff(con$res$peaks[[outer]], unique(unname(unlist(con$res$peaks[inner]))))
      excl.peaks <- con$res$peaks[[outer]]
      
      # create self-intersecting boundary for plotting or gating
      excl.bound <- innerPolyExclude(con$res$boundary[outer], con$res$boundary[inner])
      excl.names <- paste(c("*", excl.peaks), collapse = "\u00B7")
      
      # get index of exclusion if already present
      excl.idx <- which(rownames(con$res) %in% excl.names)
      if (!length(excl.idx)) {
        # add to res list
        con$res <- rbind(con$res, do.call("rbind", stats::setNames(list(tmp = list(peaks = excl.peaks, 
          p.value = NA_real_, boundary = excl.bound, bin.idx = excl.bins, 
          num.bins = length(excl.bins), eff.bins = NA_real_)), excl.names)))
        excl.idx <- nrow(con$res)
        
        # and store indices according to res list
        sel.comb[[names(con$res$bin.idx[outer])]] <- list(outer = outer, 
          inner = c(inner, excl.idx))
      }
    } else {
      sel.comb[[names(con$res$bin.idx[c])]] <- list(outer = outer, inner = inner)
    }
  }
  
  con$comb$sel <- as.data.frame(do.call("rbind", sel.comb))
  # just in case
  colnames(data) <- gsub("_|-|\\+|&|/|\\|| \\(v\\)", "", colnames(data))
  
  con$comb$res <- lapply(con$comb$sel$inner, function(inner) {
    combCompare(res = con$res, data = data, sel = list(inner), delta = delta, 
      beta = beta, threshold = threshold)
  })
  
  con$comb$data <- data
  
  return(con)
  
}
