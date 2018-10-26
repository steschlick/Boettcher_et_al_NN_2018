#' @title Assign Phenotypes
#' @description categorizes selected subsets according to 1D 'cutoffs' and assigns phenotype labels
#' @param con object returned by \code{\link{contourMerge}} or \code{\link{findCutoffs}}
#' @param cutoffs list of vectors with cutoff-values per marker
#' @param subsets indices of subsets to select
#' @param beta numeric probability value in [0,1] for symmetric quantiles
#' @param verbose flag to print out phenotypes, Default: TRUE
#' @return A list with fields:
#' \describe{
#'   \item{\code{phenocode}}{matrix with char encoded phenotypes}
#'   \item{\code{phenotype}}{phenotypes comprising all markers}
#'   \item{\code{nonredundant}}{phenotypes excluding markers that show the same phenotype among all subsets}
#' }
#' @seealso 
#'  \code{\link{phenoCutoff}}
#' @rdname cuts2Phenotype
#' @export 
#' @importFrom knitr kable
#' @importFrom stats setNames
cuts2Phenotype <- function(con, cutoffs, subsets, beta, verbose = TRUE) {
  
  if (missing(beta)) {
    beta <- con$beta
  }
  if (!all(names(cutoffs) %in% colnames(con$comb$data))) {
    stop("cutoff names don't match parameter names")
  }
  
  if (missing(subsets)) {
    subsets <- con$merge
    names(subsets) <- rownames(con$res)[subsets]
  }
  
  data <- con$comb$data[, names(cutoffs)]
  phenocode <- t(sapply(subsets, function(cl) {
    sapply(names(cutoffs), function(ch) {
      phenoCutoff(x = data[con$res$bin.idx[[cl]], ch], cutoffs = cutoffs[[ch]], 
        beta = beta)
    })
  }))
  rownames(phenocode)[grepl("\\*\u00B7root", rownames(phenocode))] <- "*\u00B7root"
  
  discr.code <- apply(phenocode, 2, function(x) !all(x == x[1]) & !any(x == ""))
  nonred.code <- apply(phenocode, 2, function(x) !all(x == x[1]))
  
  full.phenotypes <- apply(apply(phenocode, 1, function(x) {
    paste(ifelse(nzchar(x), colnames(phenocode), ""), x, sep = "")
  }), 2, paste, collapse = "")
  nonred.phenotypes <- apply(apply(phenocode[, nonred.code], 1, function(x) {
    paste(ifelse(nzchar(x), colnames(phenocode[, nonred.code]), ""), x, sep = "")
  }), 2, paste, collapse = "")
  
  uniqu.code <- phenocode[, discr.code, drop = F]
  names.ident <- apply(unique(uniqu.code), 1, function(x) {
    which(apply(uniqu.code, 1, function(y) identical(x, y)))
  })
  if (is.matrix(names.ident)) {
    names.ident <- apply(names.ident, 2, function(x) paste(rownames(uniqu.code)[x], 
      collapse = ","))
  } else {
    names.ident <- names(names.ident)
  }
  uniqu.code <- unique(uniqu.code)
  rownames(uniqu.code) <- names.ident

  # TODO: use generics (show, print, summary)
  if (verbose) {
    knitr.kable.NA <- getOption("knitr.kable.NA")
    options(knitr.kable.NA = "")
    cat("\nSubset phenotypes:")
    print(knitr::kable(t(phenocode[, nonred.code, drop = F]), align = "c"))
    cat("\nwhith all being equal in:")
    print(knitr::kable(apply(phenocode[, !nonred.code, drop = F], 2, unique), align = "c"))
    cat("\nUnique & distinct phenotypes:")
    print(knitr::kable(t(uniqu.code), align = "c"))  # 
    options(knitr.kable.NA = knitr.kable.NA)
  }
  return(stats::setNames(list(phenocode, full.phenotypes, nonred.phenotypes), c("phenocode", 
    "phenotype", "nonredundant")))
}
