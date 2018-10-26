#' @title Subset-discriminating Cutoffs for Projection-based Gating
#' @description computes cutoffs between subset "quantiloids" and finds possible projections 
#'        that allow discrimination with at most 3 cutoffs per marker
#' @param con object returned by \code{\link{contourMerge}} or \code{\link{exprs2con}}
#' @param delta numeric threshold for robust effect size \eqn{\Delta}, needs to be provided if \code{recomp=TRUE}
#' @param beta numeric threshold for non-overlapping quantiles
#' @param threshold numeric for minimally required distance between two cutoffs (and between two distributions if \code{recomp=TRUE})
#' @param recomp logical flag whether to recompare subset separation using \code{\link{combCompare}} 
#'        according to newly provided values for \code{delta}, \code{beta}, and \code{threshold}, Default: FALSE
#' @param find.all logical flag whether to ignore the result of \code{\link{combCompare}}, Default: FALSE
#' @param max.cuts the maximum number of desired cutoffs per marker, for which phenotypic separation will be tested, Default: 3
#' @param verbose flag to print out phenotypes according to found cutoffs and information 
#'        about possible phenotypic separation using \code{max.cuts} cutoffs, Default: TRUE
#' @return \code{con} object with added fields:
#' \describe{
#'   \item{\code{cutoffs}}{list with cutoffs per marker element}
#'   \item{\code{subset}}{named indices of subsets that phenotypically differ from at least one other subset}
#'   \item{\code{hinges}}{list with subsets' symmetric quantiles per marker element}
#'   \item{\code{mat.ind}}{matrix with indices of compared subsets in two columns}
#'   \item{\code{keep.mat}}{logical matrix indicating subset separation with markers in columns and subset comparisons in rows}
#'   \item{\code{delta.mat}}{matrix of effect sizes with markers in columns and subset comparisons in rows}
#'   \item{\code{sep.mat}}{square matrix summarizing which subsets are distinguishable}
#'   \item{\code{init}}{matrix with subset phenotypes according to all found cutoffs}
#'   \item{\code{reduc}}{\describe{
#'                  \item{\code{max.cuts}}{maximum number of desired cutoffs per marker}
#'                  \item{\code{cutoff.comb}}{list of per-marker cutoff combinations}
#'                  \item{\code{cutoff.comb.ind}}{indices mapping cutoff combinations to marker columns}
#'                  \item{\code{cutoff.comb.keep}}{constraints matrix for subset comparisons}
#'                  \item{\code{cutoff.comb.constraint}}{constraints matrix for cutoff combinations}
#'                  \item{\code{cutoff.comb.phenocode}}{phenotype matrix for each subset in rows and cutoff combination in columns}
#'                }}
#' }
#' @details Only markers that provide separation of at least two subsets according to the \code{delta} and \code{beta}  
#'          limits as defined in \code{\link{contourMerge}} or \code{\link{exprs2con}} are being considered. This can  
#'          be overridden or controlled by \code{recomp=TRUE} and providing new threshold values for \code{delta} and \code{beta}.
#'          Detection of cutoffs depends on non-overlapping quantiles only and can therefore be tuned 
#'          by changing the value for \code{beta}. Setting \code{find.all=TRUE} will compute cutoffs
#'          inbetween two subsets even if they did not separate according to the \code{delta} and \code{beta}
#'          limits during comparison. Depending on the number and distribution of subsets, there might be an infeasibly large number
#'          of cutoffs per marker detected, and on the other hand, considerable redundancy in the phenotype definitions. The function
#'          performs an initial check whether subsets are discriminable with \code{max.cuts} or at most 3 cutoffs, a number that 
#'          usually suffices for subset phenotyping (see \code{\link{phenoQuant}}). A final set of markers and cutoffs can then   
#'          be computed using \code{\link{phenoReduce}}. 
#' @examples 
#' \dontrun{
#' set.seed(123) 
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2, 
#'                   dimnames=list(NULL, c("A", "B"))),
#'            matrix(rnorm(36, mean = 1, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = -1, sd = 0.3), ncol = 2))
#' BW <- chooseBandWidth(x) 
#' 
#' # get density estimate
#' fhat <- BW$fhats[[ BW$Hindex[1] ]]
#' 
#' # get contour list
#' cs <- contourLoc(fhat=fhat, offset = 0.01, n.levels = 20, peaks2merge = 3, minBins = 4,
#'                  keep.level.by = "max.bin.num", plotit = TRUE)
#' 
#' # remove contours containing too few observations
#' css <- shrinkSpur(cs) 
#' 
#' # compare and merge contour-gated subsets
#' csm <- contourMerge(css, keep.outer=F, keep.inner=F, data=x, 
#'                     delta=2, beta=0.15, threshold=0.2)
#' plotContours(csm, col="grey", pch="+", cex=0.5, lwd=2.5, fhat.levels=T, labels="peak.ids")
#' 
#' # compute cutoffs
#' cscut <- findCutoffs(csm, verbose=FALSE)
#' plotPheno(cscut, comp.sel=cscut$subset, plot.cutoffs=TRUE)
#' 
#' # reduce number of markers and cutoffs 
#' csred <- phenoReduce(cscut, target.num.marker = NA, b = 0.5)
#' plotPheno(csred, comp.sel=csred$subset, plot.cutoffs=TRUE)
#' 
#' # require 1 cutoff per marker
#' cscut <- findCutoffs(csm, max.cuts=1, verbose=FALSE)
#' csred <- phenoReduce(cscut, target.num.marker = NA, b = 0.5)
#' plotPheno(csred, comp.sel=csred$subset, plot.cutoffs=TRUE)
#' }
#' @seealso 
#'  \code{\link{chooseBandWidth}}
#'  \code{\link{contourLoc}}
#'  \code{\link{shrinkSpur}}
#'  \code{\link{contourMerge}}
#'  \code{\link{phenoReduce}}
#'  \code{\link{combCompare}}
#'  \code{\link{cuts2Phenotype}}
#'  \code{\link{plotPheno}}
#'  \code{\link{plotContours}}
#'  \code{\link{phenoQuant}}
#' @rdname findCutoffs
#' @export 
#' @importFrom stats quantile setNames
#' @importFrom utils combn
#' @importFrom knitr kable
findCutoffs <- function(con, delta, beta, threshold, recomp = FALSE, 
                        find.all = FALSE, max.cuts = 3, verbose = TRUE) {
  
  if (recomp) {
    con$comb$res <- lapply(con$comb$sel$inner, function(inner) {
      combCompare(res = con$res, data = con$comb$data, sel = list(inner), delta = delta, 
        beta = beta, threshold = threshold)
    })
    con$beta <- beta
    con$delta <- delta
    con$threshold <- threshold
  }
  
  if (length(con$comb$res) > 1) {
    comb.idx <- unlist(sapply(con$comb$res, "[[", 1), recursive = FALSE)
    keep.mat <- do.call("rbind", sapply(con$comb$res, "[[", 2))
    delta.mat <- do.call("rbind", sapply(con$comb$res, "[[", 3))
    beta.mat <- do.call("rbind", sapply(con$comb$res, "[[", 4))
  } else {
    comb.idx <- con$comb$res[[1]]$idx
    keep.mat <- con$comb$res[[1]]$keep
    delta.mat <- con$comb$res[[1]]$delta
    beta.mat <- con$comb$res[[1]]$beta
  }
  
  if (!missing(delta) && !recomp) {
    keep.mat <- abs(delta.mat) >= delta
    con$delta <- delta
  }
  if (missing(beta)) {
    beta <- con$beta
  } else {
    con$beta <- beta
  }
  if (missing(threshold)) {
    threshold <- con$threshold
  } else {
    con$threshold <- threshold
  }
  
  # use only relevant markers
  marker.keep <- apply(keep.mat, 2, any)
  if (!any(marker.keep)) {
    stop("at least two separable subsets are needed, \nconsider to decrease beta and delta to find more cutoffs \nor include more markers by setting recomp=TRUE")  	
  }
  # use only relevant subset comparisons
  comb.keep <- apply(keep.mat, 1, any)
    
  # update
  data <- con$comb$data[, marker.keep]
  keep.mat <- keep.mat[comb.keep, marker.keep]
  delta.mat <- delta.mat[comb.keep, marker.keep]
  
  comb.subset <- sort(unique(unlist(comb.idx[comb.keep])))
  names(comb.subset) <- rownames(con$res)[comb.subset]
  
  mat.ind <- do.call("rbind", comb.idx[comb.keep])
  # assign indices of used subsets
  mat.ind[, 1] <- match(mat.ind[, 1], comb.subset)
  mat.ind[, 2] <- match(mat.ind[, 2], comb.subset)
  # setup lookup matrix for marker separation (based on beta & delta) between each
  # pairwise subsets as precalculated by combCompare()
  keeps <- keepm2list(keep.mat, mat.ind, names(comb.subset))
  
  # TODO: pass length.out as arg
  hinges <- sapply(colnames(data), function(m) {
    aperm(sapply(con$res$bin.idx[comb.subset], function(b) {
      cbind(lower = stats::quantile(data[b, m], probs = seq(beta, 0, length.out = 10)), 
        upper = stats::quantile(data[b, m], probs = seq(c(1 - beta), 1, length.out = 10)))
    }, simplify = "array"), c(3, 2, 1))
  }, simplify = FALSE)
  
  ## compute the cutoffs
  cutoffs <- mapply(FUN = hingeQuantileCuts, hings = hinges, keep = keeps, threshold = threshold, 
    find.all = find.all, SIMPLIFY = FALSE)
  # reduce to beta-quantiles
  hinges <- lapply(hinges, function(h) h[, , 1])
  
  con$update <- unname(comb.subset[!grepl("\\*", names(comb.subset))])
  con$keep.mat <- keep.mat
  con$delta.mat <- delta.mat
  
  # initial phenotypes, i.e. based on all found cutoffs
  phenoList <- cuts2Phenotype(con, cutoffs, comb.subset, verbose = FALSE)

  ## setup for phenoReduce ##
  ###########################

  # setup combinations of markers with limited number of cutoffs to be used in ILP,
  # create constraints as well
  max.cuts <- min(max.cuts, 3)
  data <- data[, colnames(keep.mat)]
  
  # create per-marker cutoff combinations if number of cutoffs exceeds the limit
  # (max.cuts, will always be no more than 3, to allow managable manual
  # gating/visualization)
  cutoff.comb <- list()  # holds cutoff combos
  
  for (num.cuts in 1:max.cuts) {
    cutoff.comb <- c(cutoff.comb, unlist(lapply(cutoffs, function(x) {
      l <- length(x)
      if (l > num.cuts) {
        nms <- utils::combn(l, num.cuts, simplify = FALSE)
        # names indicate which cutoffs were combined
        nms <- lapply(nms, function(n) paste0("_", paste(n, collapse = "\u00B7")))
        stats::setNames(utils::combn(x, num.cuts, simplify = FALSE), nms)
      } else list(x)
    }), recursive = FALSE))
  }
  
  # some formatting
  cutoff.comb <- cutoff.comb[!duplicated(cutoff.comb)]
  cutoff.comb <- cutoff.comb[!sapply(cutoff.comb, function(x) any(is.na(x)))]
  cutoff.comb <- cutoff.comb[order(names(cutoff.comb))]
  cutoff.comb.names <- strsplit(names(cutoff.comb), "._", fixed = TRUE)
  cutoff.comb.names <- sapply(cutoff.comb.names, "[[", 1)
  cutoff.comb.ind <- match(cutoff.comb.names, names(cutoffs))
  
  ## projecting the constraints ##
  ################################
  
  # get per-marker subset identity matrix, i.e. whether each pair of subsets are
  # distinguishable based on cutoff-phenotypes which are represented as logical
  # vectors through phenogate()
  cutoff.comb.ident <- mapply(FUN = function(hings, cuts) {
    r <- range(c(hings, cuts))
    ph <- apply(hings, 1, function(h) phenogate(q = h, c = cuts, r = r))
    !apply(ph, 2, function(p) apply(sweep(ph, 1, p, "&"), 2, any))
  }, hings = hinges[cutoff.comb.ind], cuts = cutoff.comb, SIMPLIFY = FALSE)
  # convert to long format
  cutoff.comb.keep <- keeps2mat(cutoff.comb.ident, mat.ind)
  rownames(cutoff.comb.keep) <- rownames(mat.ind)
  
  # update to keep only feasible comparisons
  comb.keep <- apply(cutoff.comb.keep, 1, any)
  
  if (sum(comb.keep) < 2) {
    stop("at least two separable subsets are needed, \nconsider to decrease beta and delta to find more cutoffs \nor include more markers by setting recomp=TRUE")
  }
  sep.mat <- keepm2list(matrix(comb.keep, dimnames = list(names(comb.keep), c("separ"))), 
    mat.ind, names(comb.subset))[[1]] * 1
  diag(sep.mat) <- NA
  dimnames(sep.mat) <- list(names(phenoList$phenotype), names(phenoList$phenotype))
  
  # TODO: use generics (show, print, summary)
  if (verbose) {
    knitr.kable.NA <- getOption("knitr.kable.NA")
    options(knitr.kable.NA = "")
    cat("\nPhenotypes according to computed cutoffs:")
    print(knitr::kable(data.frame(phenotype = phenoList$phenotype)))
    cat(paste0("\nusing ", length(colnames(keep.mat)), 
      " relevant markers:\n"))
    cat(paste(colnames(keep.mat), collapse = ", "))
    cat(paste0("\n\nPairwise separation using ", max.cuts," cutoff(s) per marker:"))
    print(knitr::kable(data.frame(distinction = ifelse(comb.keep, "yes", "no"))))
    cat("\nPhenotypic separation matrix:")
    print(knitr::kable(sep.mat))  # , align='c'
    options(knitr.kable.NA = knitr.kable.NA)
    if (sum(sep.mat[upper.tri(sep.mat)]) < length(upper.tri(sep.mat))) {
      message("\nno distinction between all of the subsets - consider \nto decrease beta and delta to find more cutoffs \nor include more markers by setting recomp=TRUE\n")
    }
  }
  
  # update
  cutoff.comb.keep <- cutoff.comb.keep[comb.keep, ]
  keep.mat <- con$keep.mat[comb.keep, ]
  delta.mat <- con$delta.mat[comb.keep, ]
  
  # create constraints for multiple cutoff combinations if no per-marker reduction
  # is needed, this is still useful for finding reduced sets of markers
  cutoff.comb.constraint <- t(1 * sapply(names(cutoffs), function(n) cutoff.comb.names %in% 
    n, simplify = TRUE))
  atleast <- rep(1, ncol(cutoff.comb.constraint))
  cutoff.comb.constraint <- rbind(atleast, cutoff.comb.constraint)
  colnames(cutoff.comb.constraint) <- names(cutoff.comb)
  
  # phenotype for each marker, combinations if present, and subset
  cutoff.comb.phenocode <- mapply(FUN = function(hings, cuts) {
    cd <- apply(hings, 1, function(h) pheno(q = h, c = cuts))
    apply(matrix(cd, ncol = nrow(hings)), 2, npz.vec)
  }, hings = hinges[cutoff.comb.ind], cuts = cutoff.comb)
  # update subset names
  rownames(cutoff.comb.phenocode) <- names(comb.subset)
  
  # update object
  con$cutoffs <- cutoffs
  con$subset <- comb.subset
  con$hinges <- hinges
  con$mat.ind <- mat.ind
  con$keep.mat <- keep.mat
  con$delta.mat <- delta.mat
  con$sep.mat <- sep.mat
  con$init <- phenoList
  con$reduc <- stats::setNames(list(max.cuts, cutoff.comb, cutoff.comb.ind, cutoff.comb.keep, 
    cutoff.comb.constraint, cutoff.comb.phenocode), c("max.cuts", "cutoff.comb", 
    "cutoff.comb.ind", "cutoff.comb.keep", "cutoff.comb.constraint", "cutoff.comb.phenocode"))
  return(con)
}

## @title helper function
## @description converts matrix to list of square matrices
## @param keepm summary (\code{keep}, \code{beta}, or \code{delta}) matrix with 
##        markers in columns and subset comparisons in rows
## @param mat.ind integer matrix with indices of compared subsets in two columns
## @param subset.names character vector of names of compared subsets
## @return list of square matrices
## @details converts a column-marker summary (\code{keep}, \code{beta}, or \code{delta}) matrix 
##          as returned from \code{\link{combCompare}} into list of square matrices (element per marker) 
##          according to selected (compared) subsets
#' @keywords internal
keepm2list <- function(keepm, mat.ind, subset.names) {
  if (length(unique(mat.ind[1:length(mat.ind)])) != length(subset.names)) {
    stop("check array index")
  }
  sapply(colnames(keepm), function(m) {
    keeps.mat <- diag(FALSE, length(subset.names))
    dimnames(keeps.mat) <- list(subset.names, subset.names)
    keeps.mat[rbind(mat.ind, mat.ind[, 2:1])] <- keepm[, m]
    keeps.mat
  }, simplify = FALSE)
}

## @title helper function
## @description converts \code{keep}-square matrices list into marker-column summary matrix
## @param keeps list of square matrices
## @param mat.ind integer matrix with indices of compared subsets in two columns
## @return matrix with list elements (markers) in columns and subset comparisons in rows
## @details needed for setting up constraints matrix for \code{\link{lp}}
#' @keywords internal
keeps2mat <- function(keeps, mat.ind) {
  keepm <- sapply(keeps, function(x) x[mat.ind], simplify = TRUE)
  rownames(keepm) <- paste(rownames(keeps[[1]])[mat.ind[, 1]], rownames(keeps[[1]])[mat.ind[, 
    2]], sep = " vs ")
  attributes(keepm)$names <- NULL
  keepm
}


