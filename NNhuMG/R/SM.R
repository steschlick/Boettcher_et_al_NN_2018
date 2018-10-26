#' @title The Skillings-Mack test statistic
#' @description performs a nonparametric Friedmann-type test for incomplete block or repeated measures designs 
#' @param y matrix with variables in columns and observations in rows
#' @param groups a vector or factor object giving the group for the corresponding rows in y, Default: NULL
#' @param blocks a vector or factor object giving the blocks or subjects for the corresponding rows in y, Default: NULL
#' @param simulate.p.value If TRUE, an estimated p-value based on the Monte Carlo method is calculated. Default: FALSE
#' @param B number of permutations, Default: 10000
#' @param sim.SM flag whether simulation data should be returned, Default: FALSE
#' @return  A list containing elements:
#' \describe{
#'   \item{p.value}{vector with p-values}
#'   \item{statistic}{vector with Skillings-Mack statistics}
#'   \item{df}{vector with degrees of freedom}
#' }
#' If \code{simulate.p.value=TRUE}: a vector of approximate p-values or a list with elements: 
#' \enumerate{
#'   \item list with test results based on the chi-squared distribution
#'   \item vector of Skilings-Mack statistics
#'   \item matrix containing the permutation distribution of the Skillings-Mack statistic
#' }
#' @details The Skillings-Mack test statistic is a generalization of the statistic used in Friedman's ANOVA method and in Durbin's rank test.
#' This nonparametric statistical test is useful for the data obtained from block or repeated measures designs with missing observations  
#' occurring randomly. The p-values are based on the chi-squared distribution or estimated by Monte Carlo method. 
#' The latter is recommended for approximating p-values when there are many ties and/or small designs with missing values are conducted.
#' This function reuses the implementation in \href{https://CRAN.R-project.org/package=Skillings.Mack}{Skillings.Mack} 
#' authored by Patchanok Srisuradetchai. Code has been slightly modified to allow for simultaneous Monte Carlo permutation of multiple variables,
#' e.g. to preserve the spatial dependence across a binned histogramm (or any multivariate data structure).
#' @references J.H. Skillings and G.A. Mack (1981), On the Use of a Friedman-Type Statistic in Balanced and Unbalanced Block Designs, Technometrics 23(2) 
#' @seealso 
#'  \href{https://CRAN.R-project.org/package=Skillings.Mack}{Skillings.Mack}
#' @rdname SM
#' @export 
#' @importFrom stats pchisq
sm.test <- function(y, groups = NULL, blocks = NULL, simulate.p.value = FALSE, B = 10000, 
  sim.SM = FALSE) {
  cols <- ncol(y)
  blocks <- factor(blocks)
  k <- nlevels(blocks)
  groups <- factor(groups)
  t <- nlevels(groups)
  g <- as.numeric(groups)
  b <- as.numeric(blocks)
  d <- array(NA, dim = c(t, k, cols))
  for (i in 1:nrow(y)) {
    d[g[i], b[i], ] <- y[i, ]
  }
  if (any(is.na(groups)) || any(is.na(blocks))) {
    stop("NA's are not allowed in groups or blocks")
  }
  if (any(diff(c(nrow(y), length(groups), length(blocks))))) {
    stop("y, groups and blocks must have the same length")
  }
  for (i in 1:k) {
    if (sum(is.na(d[, i, 1])) >= (t - 1)) {
      stop("Block#", blocks[b == i], " has only one observation. Please remove this block")
    }
  }
  
  sm <- SM.test.array(d, t, k, cols)
  
  if (simulate.p.value == TRUE) {
    logic.d <- is.na(d[, , 1])  # same over all cols
    num.missing <- matrix(NA, nrow = k, ncol = 1)
    for (i in 1:k) {
      num.missing[i, 1] <- sum(logic.d[, i])
    }
    
    simulated.SM <- vector("list", length = B)
    
    # permute the raw data (not the ranks as in original implementation)
    for (b.num in 1:B) {
      d.sim <- array(NA, dim = c(t, k, cols))
      for (i in 1:k) {
        if (num.missing[i, 1] == 0) {
          rand <- sample.int(t, size = t, replace = FALSE)
          d.sim[, i, ] <- as.matrix(d[rand, i, ])
        }
        if (num.missing[i, 1] > 0) {
          rand <- sample.int(t - num.missing[i, 1], size = t - num.missing[i, 
          1], replace = FALSE)
          d.sim[!logic.d[, i], i, ] <- as.matrix(d[!logic.d[, i], i, ])[rand, 
          ]
        }
      }
      simulated.SM[[b.num]] <- SM.test.array(d.sim, t, k, cols)$T
    }
    sim.Perm <- do.call("rbind", simulated.SM)
    sim.pval <- apply(t(sim.Perm) >= sm$T, 1, mean)
  }
  
  if (simulate.p.value) {
    if (sim.SM) {
      return(list(sm, simulated.SM, sim.Perm))
    } else return(sim.pval)
  } else if (!simulate.p.value) {
    return(list(p.value = stats::pchisq(sm$T, df = sm$rank.cov, lower.tail = FALSE), 
      statistic = sm$T, df = sm$rank.cov))
  }
}

## @title helper function
## @description calculates the actual test statistic
## @param d numeric array of dim c(t, k, cols) containing missing data
## @param t number of groups
## @param k number of blocks
## @param cols number of variables
## @return A list with elements:
## \describe{
##   \item{\code{T}}{vector with statistics}
##   \item{\code{rank.cov}}{vector with degrees of freedom}
## }
#' @importFrom MASS ginv
#' @importFrom matrixcalc matrix.rank
SM.test.array <- function(d, t, k, cols) {
  # helper function
  transform.matrix <- function(mat) {
    num.of.rows <- nrow(mat)
    num.of.cols <- ncol(mat)
    output <- matrix(rep(0, num.of.rows * num.of.rows), c(num.of.rows, num.of.rows))
    for (i in seq(1, num.of.rows - 1)) {
      for (j in seq(i + 1, num.of.rows)) {
        output[i, j] <- -sum(mat[i, ] * mat[j, ])
        output[j, i] <- output[i, j]
      }
    }
    for (i in 1:num.of.rows) {
      output[i, i] <- (-1) * sum(output[, i])
    }
    output
  }
  y.rank.NA <- array(NA, dim = c(t, k, cols))
  y.rank <- array(NA, dim = c(t, k, cols))
  T.v <- numeric(cols)
  rank.cov.v <- numeric(cols)
  # loop over variables
  for (n in 1:cols) {
    for (i in 1:k) {
      y.rank.NA[, i, n] <- rank(d[, i, n], ties.method = c("average"), na.last = "keep")
      y.rank[, i, n] <- ifelse(is.na(y.rank.NA[, i, n]) == TRUE, (sum(!is.na(y.rank.NA[, 
        i, n])) + 1)/2, y.rank.NA[, i, n])
    }
    A <- matrix(NA, nrow = 1, ncol = t)
    for (j in 1:t) {
      Ai <- matrix(NA, nrow = 1, ncol = k)
      for (i in 1:k) {
        Ai[1, i] <- sqrt(12/(sum(!is.na(y.rank.NA[, i, n])) + 1)) * (y.rank[j, 
          i, n] - (sum(!is.na(y.rank.NA[, i, n])) + 1)/2)
      }
      A[1, j] <- sum(Ai)
    }
    CovMat <- transform.matrix(!is.na(y.rank.NA[, , n]))
    Cov.Inv <- MASS::ginv(CovMat)
    rank.cov.v[n] <- matrixcalc::matrix.rank(Cov.Inv)
    T.v[n] <- A %*% Cov.Inv %*% t(A)
  }
  return(list(T = T.v, rank.cov = rank.cov.v))
}
