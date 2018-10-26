#' @title Quantile Gating
#' @description Internal function that calculates cutoffs between lower and upper hinge of any two subsets if possible
#' @param hings array of dim c(k, 2, 10) with precomputed hinges for a given marker, 
#'        i.e. a sequence of 10 upper and lower quantiles for each of k subsets
#' @param keep logical matrix indicating which combination of two subsets 
#'        differ in expression of what marker based on separation thresholds for \code{beta} 
#'        and \code{delta} as a result of \code{\link{combCompare}}
#' @param threshold numeric for minimally required distance between two cutoffs
#' @param find.all logical flag whether to ignore \code{keep}, Default: FALSE
#' @return numeric vector of cutoff-values
#' @details This function is used in \code{\link{findCutoffs}}.
#' @rdname hingeQuantileCuts
hingeQuantileCuts <- function(hings, keep, threshold, find.all = FALSE) {
  
  di <- dim(hings[, 1, ])  # subsets * quantiles
  hod <- apply(hings, 3, function(h) outer(h[, 1], h[, 2], "-"))  # lower - upper
  dim(hod) <- c(di[1], di)
  hom <- apply(hings, 3, function(h) outer(h[, 1], h[, 2], "+")/2)  # the mean
  dim(hom) <- c(di[1], di)
  hom[hod < 0] <- Inf  # no negatives and take ...
  hod[hod < 0] <- Inf
  if (!find.all) {
    hom[!keep] <- Inf  # ...separation threshold into account
    hod[!keep] <- Inf
  }
  ci <- apply(hod[, , 1], 1, which.min)  # closest hinges pair index
  up.lo <- cbind(seq_along(ci), ci)  # keep ids of hinge pairs
  qi <- diag(apply(hod[, ci, ], c(1, 2), which.min))  # index of mean
  up.lo <- cbind(up.lo, qi)
  ct <- hom[up.lo]  # cutoffs as mean between lower and closest upper
  rmv <- is.finite(ct)  # remove negatives and those below separation threshold
  if (!any(rmv)) {
    return(NA_real_)
  } else {
    up.lo <- up.lo[rmv, , drop = FALSE]
    ct <- ct[rmv]
    # collapse redundant cuts, i.e. ...
    ccm <- sapply(1:length(ct), function(i) {
      # ...those inbetween common hinge pairs at beta
      !pheno(c(hings[up.lo[i, 2], , 1], hings[up.lo[i, 1], , 1]), ct)
    })
    if (is.matrix(ccm)) {
      diag(ccm) <- TRUE  # just in case
      # if all cuts are inbetween hinge pairs take mean, if any is not, take the
      # other's mean
      ct <- unique(sapply(1:length(ct), function(i) {
        mean(ct[apply(ccm[, t(ccm)[i, ], drop = F], 1, all)])
      }))
    }
    # finally average cutoffs which are too close
    ct <- sort(unique(ct))
    ctd <- diff(ct)
    ctm <- which(ctd < threshold)
    while (length(ctm)) {
      i <- which.min(ctm)
      ct[c(i, i + 1)] <- mean(ct[c(i, i + 1)])
      ct <- sort(unique(ct))
      ctd <- diff(ct)
      ctm <- which(ctd < threshold)
    }
    return(ct)
  }
}



