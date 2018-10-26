#' @title Expression Data to Contours
#' @description contour-gating of raw expression data or phenotyping of clustering results
#' @param con object returned by \code{\link{contourMerge}} or \code{\link{findCutoffs}}
#' @param data expression matrix
#' @param markers to be used, vector of names matching columns in \code{data}  
#' @param clus.res vector containing class assigments, of the same length as number of rows in \code{data}
#' @param delta numeric threshold for robust effect size \eqn{\Delta}
#' @param beta numeric threshold for non-overlapping quantiles
#' @param threshold numeric for distance (at asinh scale) between two distributions to avoid unintentional 
#'        large \eqn{\Delta} values when expression levels are near zero
#' @param dwnsmpl number of target events to be sampled from each cluster or contour gate, Default: 512
#' @return A list object with following updated fields:
#' \describe{
#'   \item{\code{res}}{\describe{
#'                        \item{\code{peaks}}{list of indices or class for each contour or cluster}
#'                        \item{\code{bin.idx}}{list of indices of observations within each contour or cluster}
#'                      }}
#'   \item{\code{update}}{indices of contours or clusters, corresponding to rows in \code{res}}
#'   \item{\code{merged}}{indices of kept contours}
#'   \item{\code{data}}{downsampled expression matrix}
#'   \item{\code{comb}}{\describe{
#'                        \item{\code{sel}}{data.frame containing indices for each set of inner and outer contours}
#'                        \item{\code{res}}{list of results for subset or cluster comparisons as returned by \code{\link{combCompare}}}
#'                      }}
#'   \item{\code{eprs.data}}{}
#' }
#' @details This function may be used to update a contour list object to use event-level (instead of bin-level) data or convert results 
#'       obtained from clustering into a 'contour' object in order to perform phenotypic comparisons of clusters.
#' @rdname exprs2con
#' @export 
#' @importFrom splancs inout
exprs2con <- function(con, data, markers, clus.res, delta, beta, threshold, dwnsmpl = 512) {
  
  # no more because of heavy base outer() use for HL which is O(n^2)
  MAXSAM <- 1024L
  dwnsmpl <- min(dwnsmpl, MAXSAM)
  
  exprs.data <- data
  
  if (missing(clus.res) && !all(c("tSNE1", "tSNE2") %in% colnames(data))) {
    stop("need tSNE channels")
  } else {
    
    if (missing(con) && !missing(clus.res)) {
      # setup 'contour' object
      con <- list()
      subs <- sort(unique(clus.res))
      subs.name <- paste("cl_", subs, sep = "")
      peaks <- as.list(subs)
      tsne.subs <- sapply(subs, function(cl) clus.res %in% cl)
      tsne.subs.counts <- apply(tsne.subs, 2, sum)
      bin.idx <- apply(tsne.subs, 2, which)
      if (is.matrix(bin.idx)) {
        bin.idx <- split(bin.idx, rep(1:ncol(bin.idx), each = nrow(bin.idx)))
      }
      names(peaks) <- names(bin.idx) <- subs.name
      res <- data.frame(peaks = rep(NA, length(subs)), bin.idx = rep(NA, length(subs)), 
        row.names = subs.name)
      res$peaks <- peaks
      res$bin.idx <- bin.idx
      con$res <- res
      con$update <- con$merge <- seq_len(nrow(res))
      con$comb <- list()
      subs <- seq_along(subs)
      sel <- data.frame(outer = NA, inner = NA, row.names = "root")
      sel$outer <- list(root = NA)
      sel$inner <- list(root = subs)
      con$comb$sel <- sel
      
    } else {
      # contour gates
      cons <- !grepl("\\*\u00B7root", rownames(con$res))
      subs <- c(which(cons), which(!cons))
      tsne.subs <- sapply(con$res$boundary[cons], function(poly) splancs::inout(exprs.data[, 
        c("tSNE1", "tSNE2")], poly))
      root <- !apply(tsne.subs, 1, any)
      tsne.subs <- cbind(tsne.subs, root)
      tsne.subs.counts <- apply(tsne.subs, 2, sum)
    }
  }
  
  down.sample <- matrix(FALSE, nrow = nrow(data), ncol = ncol(tsne.subs))
  for (i in seq_along(tsne.subs.counts)) {
    if (tsne.subs.counts[i] > dwnsmpl) {
      ran <- sample.int(tsne.subs.counts[i], dwnsmpl)
      down.sample[tsne.subs[, i], i][ran] <- TRUE
    } else {
      down.sample[, i] <- tsne.subs[, i]
    }
  }
  keep <- apply(down.sample, 1, any)
  data <- data[keep, markers]
  sampled.subs <- down.sample[keep, ]
  
  for (s in seq_along(subs)) {
    con$res$bin.idx[[subs[s]]] <- which(sampled.subs[, s])
  }
  
  colnames(data) <- gsub("_|-|\\+|&|/|\\|| \\(v\\)", "", colnames(data))
  if (!(missing(delta) & missing(beta))) {
    con$comb$res <- lapply(con$comb$sel$inner, function(inner) {
      combCompare(res = con$res, data = data, sel = list(inner), delta = delta, 
        beta = beta, threshold = threshold)
    })
  }
  con$comb$data <- data
  con$exprs.data <- exprs.data
  
  return(con)
  
}
