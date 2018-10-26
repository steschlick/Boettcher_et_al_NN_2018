#' @title Quantile Phenotyping
#' @description Internal functions for comparing marker expression levels between two subsets 
#'              and phenotyping based on symmetric quantiles and cutoffs
#' @param x numeric vector of expression values
#' @param y numeric vector of expression values
#' @param beta numeric probability value in [0,1], threshold for non-overlapping quantiles
#' @param th numeric threshold for distance between two distributions, applied at the same scale 
#' of marker expression
#' @param cutoffs numeric vector of cutoff-values
#' @details \code{phenoQuant} returns a sign indicating non-overlap and direction of two marker distributions.
#'    \code{phenoCutoff} encodes a subset as a string according to the relative position  
#'    of its "quantiloid" to given cutoffs (1d gates). Using common syntax, i.e. by 
#'    combinations of "-", "+", we're able to classify 2 phenotypes per marker (1 cutoff), 
#'    5 phenotypes per marker (2 cutoffs), or 9 phenotypes per marker (3 cutoffs are rarely used in cytometric immunophenotyping):
#' \enumerate{
#'   \item \describe{
#'           \item{negative}{-}
#'           \item{positive}{+}
#'         }
#'   \item \describe{
#'           \item{negative}{- -  (below 1st cutoff)}
#'           \item{low}{-  (below 2nd, spanning 1st)}
#'           \item{positive}{+  (above 1st, spanning 2nd)}
#'           \item{dim}{+ -  (between 1st and 2nd)}
#'           \item{high}{+ +  (above 2nd)}
#'         }
#'   \item \describe{
#'           \item{negative}{- - -  (below 1st)}
#'           \item{very low}{- -  (below 2nd, spanning 1st)}
#'           \item{low}{-  (below 3rd, spanning 1st and 2nd)}
#'           \item{positive}{+  (above 1st, spanning 2nd and 3rd)}
#'           \item{lowdim}{+ - -  (between 1st and 2nd)}
#'           \item{dim}{+ -  (between 1st and 3rd)}
#'           \item{highdim}{+ + -  (between 2nd and 3rd)}
#'           \item{high}{+ +  (above 2nd, spanning 3rd)}
#'           \item{very high}{+ + +  (above 3rd)}
#'         }
#'  }
#' @rdname phenoType
#' @importFrom stats quantile median
#' @keywords internal
phenoQuant <- function(x, y, beta, th) {
  qntls.x <- stats::quantile(x, probs = c(beta, 1 - beta))
  qntls.y <- stats::quantile(y, probs = c(beta, 1 - beta))
  out <- outer(x, y, "-")
  m <- stats::median(out)  # HL
  ph.cd <- sign(apply(sign(outer(qntls.x, qntls.y, "-")), 2, sum)/2)
  if ((abs(sum(ph.cd)) == 2) && (abs(m) >= th)) {
    return(ph.cd[1])
  } else return(0)
}

## @title helper function
## @description subset phenotype for a given marker based on symmetric quantiles and cutoffs
#' @rdname phenoType
#' @importFrom stats quantile
#' @keywords internal
phenoCutoff <- function(x, cutoffs, beta) {
  qntls <- stats::quantile(x, probs = c(beta, 1 - beta))
  cd <- pheno(qntls, cutoffs)
  npz.vec(cd)
}

## @title helper function
## @description sign-encoding of a subset's phenotype based on quantiles and cutoffs
## @param q numeric vector of length 2, the (p, 1-p) symmetric quantiles of a subset
## @param c numeric vector, the cutoffs
## @return numeric vector of signs indicating relative position of subset to cutoffs
#' @keywords internal
pheno <- function(q, c) {
  sign(apply(sign(outer(q, c, "-")), 2, sum)/2)
}

## @title helper function
## @description encodes the phenotype of a subset
## @param q numeric vector of length 2, the (p, 1-p) symmetric quantiles of a subset
## @param c numeric vector, the cutoffs
## @param r numeric vector of length 2, the range of all 
##          subsets' (p, 1-p) symmetric quantiles and cutoffs
## @return logical vector indicating relative position of subset to cutoffs
#' @keywords internal
phenogate <- function(q, c, r) {
  ph <- sign(apply(sign(outer(q, c(r[1] - 1, c, r[2] + 1), "-")), 2, sum)/2)
  ch <- sign(diff(-1 * ph)/2)
  cummax(ch) & rev(cummax(rev(ch)))
}

## @title helper function
## @description converts sign-encoded phenotype to char-encoded phenotype
## @param x vector of signs including NA
## @return character string
#' @keywords internal
npz.vec <- function(x) {
  paste(sapply(x, npz), collapse = "")
}

## @title helper function
## @description converts vector of signs into chars of "-", "+", or ""
## @param x vector of signs including NA
## @return character vector 
#' @keywords internal
npz <- function(x) {
  c("-", "+", "", "")[c(-1, 1, 0, NA) %in% x]
}

