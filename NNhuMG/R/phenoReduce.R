#' @title Minimal Marker Cutoff Combinations
#' @description finds low-dimensional phenotypic representations of target subsets
#' @param con object returned by \code{\link{findCutoffs}}
#' @param target.num.marker desired least number of markers, Default: NA
#' @param b weighs recall and precision in F-score, Default: 0.5
#' @return \code{con} object with updated and added fields:
#' \describe{
#'   \item{\code{cutoffs}}{list with updated cutoffs per marker element}
#'   \item{\code{targets}}{named indices of subsets that phenotypically differ from at least one other subset}
#'   \item{\code{reduc}}{\describe{
#'                  \item{\code{keep}}{logical vector indicating final marker-cutoff-combinations}
#'                  \item{\code{solutions}}{matrix containing all feasible solutions returned by \code{\link{lp}}}
#'                  \item{\code{fscore}}{matrix with computed fscores for all feasible solutions}
#'                }}
#' }
#' @details The function solves an Integer Linear Programm (ILP) to find a set of at most 3 cutoffs and a minimal and/or required 
#'          number of markers that still allow to discriminate between subsets. In the likely case that multiple solutions are found,
#'          these are used to gate the target phenotypes and compute F-scores. The optimal solution is selected based on the best average
#'          accuracy among the subsets. If target.num.marker is specified, the ILP maximizes the sum of effective differences \eqn{Delta}
#'          between the subsets and a minimal set as well as a solution with at least the required number of markers is returned. See  
#'          \code{\link{findCutoffs}} for examples. Note that currently no print, show, or summary functions are implemented.
#' @seealso 
#'  \code{\link{findCutoffs}}
#'  \code{\link[lpSolve]{lp}}
#' @rdname phenoReduce
#' @export 
#' @importFrom lpSolve lp
#' @importFrom utils head
phenoReduce <- function(con, target.num.marker = NA, b = 0.5) {
  
  required <- !is.na(target.num.marker)
  
  # need to copy
  cutoff.comb <- con$reduc$cutoff.comb
  cutoff.comb.ind <- con$reduc$cutoff.comb.ind
  cutoff.comb.keep <- con$reduc$cutoff.comb.keep
  cutoff.comb.constraint <- con$reduc$cutoff.comb.constraint
  cutoff.comb.phenocode <- con$reduc$cutoff.comb.phenocode
  comb.subset <- con$subset
  
  # solve for minimal sets of markers, use negated ranks to select according to
  # separation
  keep.ord <- t(apply(-abs(con$delta.mat), 1, rank, ties.method = "min"))
  dimnames(keep.ord) <- dimnames(con$keep.mat)
  keep.comb.ord <- keep.ord[, cutoff.comb.ind]
  # TODO: if findCutoffs(find.all==FALSE) make sure that 
  # all(apply(!!keep.solve, 1, any)) 
  # keep.comb.mat <- 1 * (con$keep.mat[ , cutoff.comb.ind] & cutoff.comb.keep)
  keep.comb.mat <- 1 * cutoff.comb.keep
  
  if (required) {
    # solve for required number of markers
    target.num.marker <- min(target.num.marker, ncol(con$delta.mat))
    
    cat(paste0("\ncomputing optimal combination for ", target.num.marker, " markers\n"))
    delta.sums <- abs(con$delta.mat)[, cutoff.comb.ind]
    delta.sums[!cutoff.comb.keep] <- 0
    colnames(delta.sums) <- colnames(cutoff.comb.constraint)
    delta.solve <- 1 * cutoff.comb.keep
    delta.num.cuts <- sapply(cutoff.comb, length)
    num.cuts <- 1:con$reduc$max.cuts
    delta.sol <- vector("list", length(num.cuts))
    
    for (n in num.cuts) {
      d.con <- rbind(cutoff.comb.constraint, delta.solve)
      d.con[, delta.num.cuts > n] <- 0
      d.obj <- 1 + max(colSums(delta.sums)) - colSums(delta.sums)
      # set target minimum number of markers
      d.rhs <- c(target.num.marker, rep(1, nrow(d.con) - 1))
      m.dir <- rep("<=", nrow(cutoff.comb.constraint) - 1)
      d.dir <- c(">=", m.dir, rep(">=", nrow(delta.solve)))
      print(d.sol <- lpSolve::lp("min", d.obj, d.con, d.dir, d.rhs, all.bin = TRUE, 
        num.bin.solns = 100))
      delta.sol[[n]] <- matrix(utils::head(d.sol$solution, ncol(d.con) * d.sol$num.bin.solns), 
        nrow = d.sol$num.bin.solns, byrow = TRUE)
    }
    
    delta.sol <- do.call("rbind", delta.sol)
    keep.comb.mat[, !apply(!(!delta.sol), 2, any)] <- 0
  }
  
  solutions.l <- vector("list", max(keep.ord))
  cat("\ncalculating minimal marker combinations, starting with best separating markers\n")
  for (i in seq_along(solutions.l)) {
    keep.solve <- keep.comb.mat
    keep.solve[keep.comb.ord > i] <- 0
    
    f.con <- rbind(cutoff.comb.constraint, keep.solve)
    f.obj <- rep(1, ncol(f.con))
    # minimize objective function, i.e. minimal combo that satisfies constraints
    f.rhs <- c(1, rep(1, nrow(f.con) - 1))
    m.dir <- rep("<=", nrow(cutoff.comb.constraint) - 1)
    f.dir <- c(">=", m.dir, rep(">=", nrow(keep.solve)))
    
    cat("\n", i, "\n")
    # Please note, we need to restrict number of feasible solution for there will be
    # combinatorial explosion for many markers included in the model 
    # encountered crashes when 100 < num.bin.solns < 200 with expected solutions >> 200
    print(sol <- lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs, all.bin = TRUE, num.bin.solns = 100))
    solutions.l[[i]] <- matrix(utils::head(sol$solution, ncol(f.con) * sol$num.bin.solns), 
      nrow = sol$num.bin.solns, byrow = TRUE)
  }
  # remove infeasible solutions
  feasible <- sapply(solutions.l, function(sol) any(!(!sol)))
  if (!any(feasible)) 
    stop("try to increase delta, max number of cutoffs, or select subsets manually")
  solutions.l <- solutions.l[feasible]
  
  # check also for max number of possible sols
  solutions.l <- solutions.l[sapply(solutions.l, function(sol) nrow(sol) < 100)]
  
  if (required) {
    solutions.l <- c(list(delta.sol), solutions.l)
  }
  
  # setup list to track (delta separation) ranks, number of solutions, number of
  # markers, additional combinations of markers/cuts will be appended, duplicates
  # removed
  sol.track <- sapply(seq_along(solutions.l), function(d.rank) {
    tmp <- cbind(d.rank, n.marker = apply(solutions.l[[d.rank]], 1, sum))
    tmp <- cbind(tmp, solution = seq_len(nrow(tmp)))
    tmp
  }, simplify = FALSE)
  
  # collapse solutions with equal number of markers (same within each)
  sol.colps <- sapply(sol.track, "[", 1, 2)
  sol.track <- lapply(unique(sol.colps), function(sc) {
    do.call("rbind", sol.track[sol.colps %in% sc])
  })
  solutions <- lapply(unique(sol.colps), function(sc) {
    do.call("rbind", solutions.l[sol.colps %in% sc])
  })
  # remove duplicates (caused by replicated solutions)
  sol.track <- mapply(function(sol.h, sols) {
    sol.h[!duplicated(sols), , drop = FALSE]
  }, sol.track, solutions, SIMPLIFY = FALSE)
  
  solutions <- lapply(solutions, function(sols) sols[!duplicated(sols), , drop = FALSE])
  
  # take resulting marker combinations and pick the one with best average f-score
  # among the subsets.
  sol.track <- do.call("rbind", sol.track)
  all.sols <- do.call("rbind", solutions)
  colnames(all.sols) <- colnames(cutoff.comb.constraint)
  
  # rename
  names(cutoff.comb) <- colnames(keep.solve)
  rownames(cutoff.comb.phenocode) <- names(comb.subset)
  
  # create matrix of true positives (sampled from raw data)
  tsne.trupos <- sapply(con$res$bin.idx[comb.subset], function(b) seq_len(nrow(con$comb$data)) %in% 
    b)
  
  colnames(tsne.trupos) <- names(comb.subset)
  # exclude the exclusion (i.e. NOT-subset-gated cells/bins)
  tsne.trupos <- tsne.trupos[, !grepl("\\*\u00B7root", names(comb.subset))]
  cutoff.comb.phenocode <- cutoff.comb.phenocode[!grepl("\\*\u00B7root", names(comb.subset)), 
    ]
  sol.phenotypes <- sapply(seq_len(nrow(all.sols)), function(l) {
    cutoff.comb.phenocode[, !(!all.sols[l, ]), drop=FALSE]
  }, simplify = FALSE)
  sol.cutoff.comb <- sapply(seq_len(nrow(all.sols)), function(l) {
    cutoff.comb[!(!all.sols[l, ])]
  }, simplify = FALSE)
  
  # compute f-measure for each gated phenotype/solution
  fms.sols <- sapply(seq_len(nrow(all.sols)), function(l) {
    phenocode <- cutoff.comb.phenocode[, !(!all.sols[l, ]), drop=FALSE]
    dat <- con$comb$data[, names(cutoff.comb[!(!all.sols[l, ])]), drop=FALSE]
    sapply(rownames(phenocode), function(s) {
      gated <- apply(sapply(colnames(dat), function(m) {
        testCutoffs(x = dat[, m], pheno = phenocode[s, m], cutoffs = cutoff.comb[!(!all.sols[l, 
          ])][[m]])
      }), 1, all)
      fmeasure(pred = gated, true = tsne.trupos[, s], b = b)
    })
  })
  
  # w/o any gating
  fms.null <- apply(tsne.trupos, 2, function(tr) {
    fmeasure(pred = !logical(nrow(con$comb$data)), true = tr, b = b)
  })
  # max per subset
  fms.max <- apply(fms.sols, 1, max)
  fms.opt <- apply(fms.sols, 1, function(x) which(x == max(x)))
  fms.opt <- unique(unlist(fms.opt))
  
  # calculate f-measure only within a combo so we need a hierarchy here
  hierarchy.sol <- apply(all.sols, 1, function(x) {
    apply(all.sols, 1, function(y) {
      all(x[!(!x)] == y[!(!x)])
    })
  })
  
  if (!is.matrix(hierarchy.sol)) {
  	hierarchy.sol <- as.matrix(hierarchy.sol) 
  }
  
  # calculate row max for each family 
  # fms.sols.n <- sweep(fms.sols, 1, fms.max, '/')
  fms.max.hierarchy <- apply(hierarchy.sol, 1, function(x) {
    apply(fms.sols[, x, drop = FALSE], 1, max)
  })
  colnames(fms.max.hierarchy) <- seq_len(ncol(fms.max.hierarchy))
  
  # pick the best average
  opt <- which.max(apply(fms.max.hierarchy, 2, mean))
  # find maxima for each subset in (possibly) different numbers of markers
  targets <- apply(fms.sols[, hierarchy.sol[opt, ], drop = FALSE], 1, which.max)
  
  # cbind(fms.null, hierarchy.max=fms.max.hierarchy[ , opt],
  # fms.subset.max=fms.max, fms.sols[, hierarchy.sol[opt, ], drop=FALSE])
  
  # the target marker-cutoff-combo
  cuts.target <- apply(!(!all.sols[hierarchy.sol[opt, ], , drop = FALSE]), 2, any)
  
  # check, this should never happen, for testthat
  if (required && !any(apply(delta.sol, 1, function(x) all(all.sols[hierarchy.sol[opt, 
    ], , drop = FALSE][1, ] == x)))) {
    stop("check the ILP model!")
  }
  
  target.hierarchy <- sol.phenotypes[hierarchy.sol[opt, ]]
  
  target.phenos <- sol.phenotypes[which(hierarchy.sol[opt, ])[targets]]
  # extract the target phenotyp for each subset (the 'list' diagonal)
  target.phenos <- sapply(seq_along(target.phenos), function(t) target.phenos[[t]][t, 
    , drop = FALSE], simplify = FALSE)
  target.phenos <- sapply(target.phenos, function(x) paste(ifelse(nzchar(x[1, ]), 
    colnames(x), ""), x[1, ], sep = ""), simplify = FALSE)
  names(target.phenos) <- paste("m\u00B7", rownames(fms.sols), sep = "")
  
  # phenotypes for required number of markers or for any more than the minimal
  full.phenos <- unlist(apply(target.hierarchy[[1]], 1, function(x) list(paste(ifelse(nzchar(x), 
    colnames(target.hierarchy[[1]]), ""), x, sep = ""))), recursive = FALSE)
  names(full.phenos) <- paste("f\u00B7", names(full.phenos), sep = "")
  target.phenos <- c(full.phenos, target.phenos)
 # target.phenos <- target.phenos[!duplicated(target.phenos)]

  con$cutoffs <- cutoff.comb[cuts.target]
  con$targets <- target.phenos
  con$reduc$keep <- cuts.target
  con$reduc$solutions <- all.sols  
  con$reduc$fscores <- fms.sols  
  return(con) 
}

## @title helper function
## @description simple 1D-gating
## @param x numeric vector of expression values
## @param pheno char-encoded phenotype of a given marker expression
## @param cutoffs numeric vector of cutoff-values
## @return logical vector
## @details gating without flowCore overhead. see also \code{\link{pheno}}
#' @keywords internal
testCutoffs <- function(x, pheno, cutoffs) {
  n <- length(cutoffs)  # cutoffs should be sorted already
  if (n < 1 || n > 3) {
    NULL
  } else {
    switch(n, {
      if (pheno == "") return(!logical(length(x)))
      if (pheno == "-") return(x <= cutoffs[1])
      if (pheno == "+") return(cutoffs[1] < x)
    }, {
      if (pheno == "") return(!logical(length(x)))
      if (pheno == "--") return(x <= cutoffs[1])
      if (pheno == "-") return(x <= cutoffs[2])
      if (pheno == "+-") return(cutoffs[1] < x & x <= cutoffs[2])
      if (pheno == "+") return(cutoffs[1] < x)
      if (pheno == "++") return(cutoffs[2] < x)
    }, {
      if (pheno == "") return(!logical(length(x)))
      if (pheno == "---") return(x <= cutoffs[1])
      if (pheno == "--") return(x <= cutoffs[2])
      if (pheno == "-") return(x <= cutoffs[3])
      if (pheno == "+--") return(cutoffs[1] < x & x <= cutoffs[2])
      if (pheno == "+-") return(cutoffs[1] < x & x <= cutoffs[3])
      if (pheno == "++-") return(cutoffs[2] < x & x <= cutoffs[3])
      if (pheno == "+") return(cutoffs[1] < x)
      if (pheno == "++") return(cutoffs[2] < x)
      if (pheno == "+++") return(cutoffs[3] < x)
    })
  }
}


## @title helper function
## @description computes the harmonic mean of precision and recall
## @param pred logical vector of predicted positives
## @param true logical vector of true positives
## @param b (beta-) parameter weighing precision vs. recall, Default: 1
## @return F-beta score
## @details b = 0.5 weighs recall lower than precision
#' @keywords internal
fmeasure <- function(pred, true, b = 1) {
  retrieved <- sum(pred)
  if (retrieved != 0) {
    precision <- sum(pred & true)/retrieved
  } else precision <- 0
  recall <- sum(pred & true)/sum(true)
  if ((recall != 0) && (precision != 0)) {
    Fm <- (1 + b^2) * precision * recall/((b^2) * precision + recall)
  } else Fm <- 0
  Fm
}
