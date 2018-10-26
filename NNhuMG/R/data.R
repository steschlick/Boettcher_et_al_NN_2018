#' study meta data.
#' @format A data frame with 36 rows and 11 variables:
#' \describe{
#'   \item{Sample}{sample name}
#'   \item{fcs}{fcs file name}
#'   \item{Pool}{pooled sample}
#'   \item{BC}{barcode number}
#'   \item{Barcode}{barcode sample name}
#'   \item{Donor}{donor id}
#'   \item{Region}{brain region, from which samples were isolated}
#'   \item{Sex}{donor gender}
#'   \item{Age}{donor age}
#'   \item{donor.col}{figure color code used for donor}
#'   \item{region.col}{figure color code used for brain region}
#' }
"study"

#' fcs data ranges.
#' @format A matrix with 56 rows and 2 columns
"bounds"

#' emd scores.
#' @format A square matrix with 36 rows and 36 columns
"emd"

#' p-values.
#' @format A vector of 512 Monte Carlo-simulated p-values, one for each bin
"pvals"

#' Tmax statistic.
#' @format A vector of 10,000 Monte Carlo-simulated maximal Skillings-Mack statistics, one for each resampling
"tmax"
