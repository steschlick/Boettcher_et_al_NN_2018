#' @title Contour-based Gating of Significant Subsets
#' @description locates regions of density-'tagged' observations in two dimensions
#' @param data frequency matrix used to test for group-level differences in \code{statsfun} with bins in columns and samples in rows 
#' @param fhat object returned by \code{\link{kdeDiff}} or \code{\link[ks]{kde}}
#' @param offset value in [0, 1] specifying the probability for which regions are detected, i.e. the lowest magnitude 
#'        at which contours are reported, Default: 0.1
#' @param n.levels number of levels at which contours are evaluated, Default: 50
#' @param alpha level of significance for testing group differences, Default: 0.05
#' @param statsfun a function that will take for each evaluated contour (aggre-) gated data and return a \emph{p}-value for 
#'        group-level inference (see example for details). Currently implemented for binned data only.
#' @param peaks2merge number of peaks that a contour may contain. This allows to report larger contours at lower magnitude 
#'        and is useful in the presence of spurious peaks, Default: 4
#' @param minBins minimum required number of observations that a contour must encompass to be reported, Default: 10
#' @param keep.level.by a character string of "max.bin.num", or "min.p.value", "max.bin.cor", indicating whether for each dens region
#'        the largest contour (i.e. at lowest magnitude for which it belongs to the same peak(s)) will be reported, the contour with 
#'        the most significant group-level statistic, or whether correlative structure of the data is taken into account, respectively. 
#'        In the latter, an effective number of bins determines which contour with \emph{p}-value <= \code{alpha} to report.
#' @param plotit flag to plot detected contours, Default: TRUE
#' @return  A list object with following fields:
#' \describe{
#'   \item{\code{res}}{\describe{
#'                       \item{\code{peaks}}{list of peak indices for each contour}
#'                       \item{\code{p.value}}{list of p.values or NA for each contour}
#'                       \item{\code{boundary}}{list of x- and y-coordinates for each contour polyon}
#'                       \item{\code{bin.idx}}{list of indices of observations within each contour}
#'                       \item{\code{num.bins}}{list of observation numbers within each contour}
#'                       \item{\code{eff.bins}}{list of effective observation numbers within each contour}
#'                     }}
#'   \item{\code{update}}{indices of contours, corresponding to rows in \code{res}}
#'   \item{\code{fhat}}{input object returned by \code{\link{kdeDiff}} or \code{\link[ks]{kde}}}
#'   \item{\code{peaks}}{peak points matrix with coordinates and level in columns}
#'   \item{\code{levels}}{levels at which countours were evaluated}
#'   \item{\code{data}}{}
#'   \item{\code{minBins}}{}
#'   \item{\code{peaks2merge}}{}
#'   \item{\code{offset}}{}
#'   \item{\code{alpha}}{}
#'   \item{\code{statsfun}}{}
#'   \item{\code{keep.level.by}}{}
#' }
#' @details The functionalities provided with \code{\link{contourLoc}} and accompanying functions (\code{\link{shrinkSpur}}, \code{\link{contourMerge}})
#'        integrate into a larger workflow and aim at facilitating exploration, interaction with, and interpretation of the cytometry data after statistical 
#'        testing in a 'high-resolution' setting. The input is expected to be an object returned by \code{\link{kdeDiff}}, where individual bin-coordinates
#'        (of a \code{\link[flowFP]{flowFPModel}} derived from a composite t-SNE map) as observations have been weighted by their \emph{p}-values prior   
#'        to kernel density estimation. The density estimate thus contains local information about the significance testing and this is exploited here to  
#'        locate regions encompassing significant bins. As such, \code{\link{contourLoc}} is \emph{per se} not designed for automated gating of the actual
#'        observations or events (although this works as well). This may more conveniently be achieved by \code{\link[flowStats]{curv2Filter}}. Note that 
#'        currently no print, show, or summary functions are implemented.
#' @examples 
#' \dontrun{
#' 
#' # mock 118 bin coordinates
#' set.seed(123)
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2,
#'                   dimnames=list(NULL, c("A", "B"))),
#'            matrix(rnorm(36, mean = 1, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = -1, sd = 0.3), ncol = 2))
#' plot(x) # 3 subsets
#'
#' # row indices of 2nd and 3rd subsets
#' sub2 <- 51:68
#' sub3 <- 69:118
#'
#' # mock frequency data for 118 bins and 20 samples (10 each group)
#' b <- jitter(matrix(1, nrow=20, ncol=118))
#' 
#' # add frequency difference for 2nd subset
#' r1 <- matrix(rep(rnorm(10, 0.01, 0.005) , 18), ncol=18, byrow=TRUE)
#' b[1:10, sub2] <- b[1:10, sub2] + r1 # more in group 1
#' b[11:20, sub2] <- b[11:20, sub2] - r1
#' 
#' # add frequency difference within 3rd subset
#' r2 <- matrix(rep(rnorm(10, 0.01, 0.001) , 50), ncol=50, byrow=TRUE)
#' shift <-  (x[sub3, 2] / x[sub3, 1]) < 1 
#' b[1:10, sub3[shift]] <- b[1:10, sub3[shift]] - r2[, shift]
#' b[11:20, sub3[shift]] <- b[11:20, sub3[shift]] + r2[, shift]
#' b[1:10, sub3[!shift]] <- b[1:10, sub3[!shift]] + r2[, !shift]
#' b[11:20, sub3[!shift]] <- b[11:20, sub3[!shift]] - r2[, !shift]
#' 
#' # frequency mock ready
#' b <- sweep(b, 1, rowSums(b), '/')
#' 
#' # two one-sided tests to see direction of difference
#' pv1 <- apply(b, 2, function(x) wilcox.test(x[1:10], x[11:20], alternative="greater")$p.value)
#' pv2 <- apply(b, 2, function(x) wilcox.test(x[1:10], x[11:20], alternative="less")$p.value)
#' 
#' # scale p-values to weights
#' w1 <- threshold2one(-log2(p.adjust(pv1, "BH")), -log2(0.05))
#' w2 <- threshold2one(-log2(p.adjust(pv2, "BH")), -log2(0.05))
#' 
#' # get density estimate for bin coordinates
#' BW <- chooseBandWidth(x) 
#' f0 <- BW$fhats[[ BW$Hindex[1] ]]
#' 
#' # get density-difference with 'significance weights'
#' f1 <- kdeDiff(fhat=f0, w1=w1)
#' f2 <- kdeDiff(fhat=f0, w1=w2)
#' 
#' # statsfun for testing contour-gated data 
#' statsfun <- function(data) wilcox.test(data[1:10,], data[11:20,])$p.value	
#' 
#' # get subsets more abundant in group 1
#' gr <- contourLoc(data=b, fhat=f1, offset = 0.1, statsfun=statsfun, n.levels = 20, peaks2merge = 2, 
#' minBins = 4, keep.level.by = "min.p.value", plotit = FALSE)
#' plotContours(gr, fhat.levels=T, labels="peak.ids", xlim=c(-2, 2), ylim=c(-2, 2))
#' 
#' # get subsets more abundant in group 2
#' le <- contourLoc(data=b, fhat=f2, offset = 0.1, statsfun=statsfun, n.levels = 20, peaks2merge = 2, 
#' minBins = 4, keep.level.by = "min.p.value", plotit = FALSE)
#' plotContours(le, fhat.levels=T, labels="peak.ids", xlim=c(-2, 2), ylim=c(-2, 2))
#' }
#' @seealso 
#'  \code{\link[grDevices]{contourLines}}
#'  \code{\link{shrinkSpur}}
#'  \code{\link{contourMerge}}
#'  \code{\link{plotContours}}
#'  \code{\link{plotPheno}}
#'  \code{\link{findCutoffs}}
#'  \code{\link{phenoReduce}}
#' @rdname contourLoc
#' @export 
#' @import stats
#' @importFrom grDevices contourLines
#' @import graphics
contourLoc <- function(data, fhat, offset = 0.1, n.levels = 50, alpha = 0.05, statsfun, peaks2merge = 4, minBins = 10, 
                       keep.level.by = c("max.bin.cor", "max.bin.num", "min.p.value"), plotit = TRUE) {
  
  stats <- !missing(statsfun)
  if (!stats) 
    statsfun <- NA
  # if ((keep.level.by %in% c('max.bin.cor','min.p.value')) & !missing(data)) {
  if (!missing(data)) {
    use.eff <- TRUE
  } else {
    use.eff <- FALSE
    data <- NA
  }
  
  minBins <- max(2, minBins)
  d <- stats::predict(fhat, x = fhat$x)
  offset.lvl <- stats::quantile(d[d != 0], prob = offset)
  
  # peak point matrix
  p.mat <- PeakFind(fhat, tol = offset.lvl)
  # TODO: check for peaks at the very same level, if so, add some noise?
  if (length(unique(p.mat[, 3])) != nrow(p.mat)) 
    stop("TODO!")
  if (is.na(p.mat[1])) 
    stop("no peaks detected")
  
  # need to have at least 1 level between each two peaks / local maxima  
  lvls <- c(sort(p.mat[, 3], decreasing = TRUE), offset.lvl, 0)
  delta.lvls <- lvls[-length(lvls)] - lvls[-1]
  hlf.lvls <- sort(c(lvls[length(lvls)], t(sapply(delta.lvls, "/", c(1:2))) + lvls[-1]), 
    decreasing = TRUE)
  pr.lvls <- pretty(c(sort(c(p.mat[, 3], hlf.lvls))), n = n.levels)
  lvls <- sort(unique(c(hlf.lvls, pr.lvls[-which(pr.lvls > max(hlf.lvls))])), decreasing = TRUE)
  lvls <- lvls[1:which(lvls == offset.lvl)]
  
  # ContourL holds contour lines for each level 
  # PeakLabels holds indices of peaks in each contour
  ContourL <- PeakLabels <- vector("list", length(lvls))
  
  # descend through countour levels
  for (i in 1:length(lvls)) {
    ContourL[[i]] <- grDevices::contourLines(fhat$eval.points[[1]], fhat$eval.points[[2]], 
      fhat$estimate, nlevels = 1, levels = lvls[i])
  }
  for (i in 1:length(ContourL)) {
    # collect peak ids for each contour of that level
    PeakLabels.tmp <- sapply(ContourL[[i]], function(CL) {
      list(inpoly(p.mat, CL))
    })
    # remove inverse contours, i.e. those not containing any peak
    keep <- as.logical(sapply(PeakLabels.tmp, length))
    PeakLabels[[i]] <- PeakLabels.tmp[keep]
    ContourL[[i]] <- ContourL[[i]][keep]
  }
  
  if (plotit) {
    xlim <- range(fhat$x[, 1])
    ylim <- range(fhat$x[, 2])
    graphics::plot(fhat$x, col = "grey", pch = "*", xlab = "tSNE1", ylab = "tSNE2", xlim = xlim, 
      ylim = ylim)
    graphics::contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, col = "grey", 
      lwd = 0.2, levels = lvls, drawlabels = FALSE, add = TRUE)
  }
  
  res <- list()
  
  count <- 0  # region counter
  
  for (m in 1:peaks2merge) {
    
    avlbl <- 1:nrow(p.mat)
    
    while (length(avlbl) > 0) {
      i <- avlbl[1]  # index of target peak
      
      # indices of contour levels that contain peak i and m-1 more peaks keep a list of
      # length(PeakLabels) to enable subsetting ContourL
      PLidx <- sapply(PeakLabels, function(p) {
        which(sapply(p, function(c) (i %in% c) && (length(c) == m)) == TRUE)
      }, simplify = FALSE)
      
      # two contours might contain the same peak at the same level, i.e. the innermost
      # contour is enclosed by a valley!
      PLencl <- which(sapply(PLidx, length) > 1)
      if (length(PLencl)) {
        # remove outer contours!
        for (x in PLencl) {
          PLidx[[x]] <- PLidx[[x]][which.min(sapply(ContourL[[x]][PLidx[[x]]], 
          function(CL) arpoly(CL)))]
        }
      }
      
      # keep contour for given peak get contour with maximum (effective) bin number or
      # minimum p.value
      PLtarget <- which(sapply(PLidx, length) == 1)
      
      res.tmp <- sapply(PLtarget, function(x) {
        boundary <- list(x = ContourL[[x]][[PLidx[[x]]]]$x, y = ContourL[[x]][[PLidx[[x]]]]$y)
        
        bin.idx <- inpoly(fhat$x, boundary)
        nb <- length(bin.idx)
        if (nb >= minBins) {
          if (!use.eff) {
          p.value <- NA
          eff.bins <- nb
          } else {
          bin.data <- matrix(data[, bin.idx], ncol = nb)
          if (stats) {
            # aggregate data of target bins
            data.tmp <- matrix(apply(bin.data, 1, sum), ncol = 1)
            # NA in eff.bins control signif level alpha or minBins
            p.value <- statsfun(data.tmp)
            eff.bins <- ifelse(p.value <= alpha, getEffectivBinCor(bin.data), 
            NA)
          } else {
            p.value <- NA
            eff.bins <- getEffectivBinCor(bin.data)
          }
          }
        } else {
          p.value <- NA
          eff.bins <- NA
        }
        list(peaks = PeakLabels[[x]][[PLidx[[x]]]], p.value = p.value, boundary = boundary, 
          bin.idx = bin.idx, num.bins = nb, eff.bins = eff.bins)
      }, simplify = FALSE)
      
      # NA in eff.bins control signif level alpha or minBins
      comp.bin.cor <- sapply(res.tmp, "[[", 6)
      
      keep.level.by <- match.arg(keep.level.by)
      compare.tmp <- switch(keep.level.by, max.bin.cor = comp.bin.cor, max.bin.num = {
        comp.bin.num <- sapply(res.tmp, "[[", 5)
        comp.bin.num[is.na(comp.bin.cor)] <- NA
        comp.bin.num
      }, min.p.value = {
        comp.p.value <- sapply(res.tmp, "[[", 2)
        comp.p.value[is.na(comp.bin.cor)] <- NA
        -log2(as.numeric(comp.p.value))
      })
      
      if (!all(is.na(compare.tmp))) {
        
        sel <- which(compare.tmp == max(compare.tmp, na.rm = TRUE))
        
        # if multiple max, get contour at lowest, i.e. last level
        sel.idx <- PLtarget[sel[length(sel)]]
        
        if (plotit) {
          graphics::polygon(ContourL[[sel.idx]][[PLidx[[sel.idx]]]]$x, ContourL[[sel.idx]][[PLidx[[sel.idx]]]]$y, 
          border = m, xlim = xlim, ylim = ylim, lwd = 2)
        }
        
        count <- count + 1
        
        res[count] <- res.tmp[sel[length(sel)]]
        
        avlbl <- avlbl[!(avlbl %in% PeakLabels[[sel.idx]][[PLidx[[sel.idx]]]])]
      } else avlbl <- avlbl[-1]
    }
  }
  names(res) <- sapply(sapply(res, "[[", 1), paste, collapse = "\u00B7")  # &middot;
  
  res <- as.data.frame(do.call("rbind", res))
  update <- seq_len(nrow(res))
  
  return(stats::setNames(list(res, update, fhat, p.mat, lvls, data, minBins, peaks2merge, 
    offset, alpha, statsfun, keep.level.by), c("res", "update", "fhat", "peaks", 
    "levels", "data", "minBins", "peaks2merge", "offset", "alpha", "statsfun", 
    "keep.level.by")))
}

## @description helper function
## @param x numeric matrix
## @param method character string of "pearson", "kendall", or "spearman" (Default) used to compute correlation matrix
## @return numeric number of effective bins
## @details This function computes an effective bin number based on the correlative structure of the data.
#' @importFrom stats cor
getEffectivBinCor <- function(x, method = "spearman") {
  corMatr <- stats::cor(x, method = method, use = "na.or.complete")
  nrow(corMatr) * mean(corMatr)
}



