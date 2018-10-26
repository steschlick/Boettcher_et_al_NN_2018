#' @title Contour Plot
#' @description high-level plotting function
#' @param con object returned by \code{\link{contourLoc}}, \code{\link{shrinkSpur}}, \code{\link{contourMerge}}, 
#'        \code{\link{findCutoffs}}, or \code{\link{phenoReduce}}
#' @param fhat.levels flag to plot contour levels, Default: FALSE
#' @param col color for observations points, Default: 'grey'
#' @param pch plotting 'character', i.e., symbol to use for observations points, Default: '+'
#' @param cex size of observations points, Default: 0.5
#' @param lwd width of contour lines, Default: 2.5
#' @param comp.sel indices of contours to plot, Default: FALSE
#' @param density the density of shading lines, in lines per inch, Default: NULL
#' @param filwd width of shading lines, Default: 1.5
#' @param labels how peaks shall be plotted, one of "none", "peak.pts", or "peak.ids"
#' @param cex.labels size of peak labels, Default: 0.5
#' @param ask passed to \code{\link[graphics]{par}}, Default: TRUE
#' @param ... further args passed to \code{\link[graphics]{plot}}
#' @rdname plotContours
#' @export 
#' @import graphics
plotContours <- function(con, fhat.levels = FALSE, col = "grey", pch = "+", cex = 0.5, 
  lwd = 2.5, comp.sel = FALSE, density = NULL, filwd = 1.5, labels = c("none", 
    "peak.pts", "peak.ids"), cex.labels = 0.5, ask = TRUE, ...) {
  
  if (comp.sel && !is.numeric(comp.sel) && !is.null(con$comb)) {
    sq <- seq_len(nrow(con$comb$sel))
    graphics::par(ask = ask)
    for (i in sq) {
      sel <- con$comb$sel$inner[[i]]
      plotContours(con, comp.sel = sel, fhat.levels = FALSE, col = col, pch = pch, 
        cex = cex, lwd = lwd, filwd = filwd, density = density, labels = labels, 
        cex.labels = cex.labels, ask = ask, ...)
      
    }
    graphics::par(ask = FALSE)
    return(invisible(NULL))
  }
  
  if (is.numeric(comp.sel)) {
    sel <- comp.sel
  } else {
    sel <- con$update
  }
  
  # temporary res object
  obj <- con$res[sel, ]
  # exclusion 'gates'
  ex <- grepl("\\*", rownames(obj))
  
  con.bound <- obj$boundary
  con.n <- length(con.bound)
  con.cols <- gg_color_hue(con.n)
  graphics::plot(con$fhat$x, col = col, pch = pch, cex = cex, lwd = lwd, ...)
  
  if (fhat.levels) {
    graphics::contour(con$fhat$eval.points[[1]], con$fhat$eval.points[[2]], con$fhat$estimate, 
      col = "grey", lwd = 0.2, levels = con$levels, drawlabels = FALSE, add = TRUE)
  }
  
  b_len <- seq_len(con.n)
  if (any(ex)) 
    b_len <- c(which(ex), b_len[-which(ex)])
  
  for (b in b_len) {
    if (ex[b]) {
      d <- density
      c <- con.cols[b]
    } else {
      d <- NULL
      c <- NA
    }
    
    graphics::polygon(con.bound[[b]]$x, con.bound[[b]]$y, border = con.cols[b], density = d, 
      col = c, fillOddEven = ex[b], lwd = ifelse(ex[b], filwd, lwd))
  }
  
  labels <- match.arg(labels)
  switch(labels, none = invisible(NULL), peak.pts = graphics::points(con$peaks[, 1:2], pch = 19, 
    cex = cex.labels), peak.ids = graphics::text(con$peaks[, 1:2], labels = seq(1:nrow(con$peaks)), 
    cex = cex.labels))
  
}

## @description helper function to emulate ggplot2 default (factor) color palette in base graphics
## @param n factor levels
## @return vector of hex colors
## @rdname plotContours
#' @importFrom grDevices hcl
gg_color_hue <- function(n) {
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
