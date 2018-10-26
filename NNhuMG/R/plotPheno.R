#' @title Plot Subset Marker Expression
#' @description high-level plotting function
#' @param con object returned by \code{\link{contourMerge}}, \code{\link{findCutoffs}}, or \code{\link{phenoReduce}} 
#' @param comp.sel indices of contour-gated subsets to plot, Default: FALSE
#' @param marker.ord order of markers, Default: NA uses alphabetical order
#' @param relevant.only flag whether to exclude noninformative markers, Default: FALSE
#' @param plot.cutoffs flag to plot horizontal cutoff-lines, Default: FALSE
#' @param delta.panel flag whether to plot barcharts for pairwise computed \eqn{Delta} values, Default: TRUE
#' @param max.delta.bars number of bars to plot per dodge, Default: 6
#' @param ask passed to \code{\link[graphics]{par}}, Default: TRUE
#' @param grobs flag whether to return a list of grobs without drawing, see \code{\link[gridExtra]{arrangeGrob}}, Default: FALSE
#' @seealso 
#'  \code{\link{plotContours}}
#' @rdname plotPheno
#' @export 
#' @importFrom graphics par
#' @importFrom grDevices boxplot.stats
#' @importFrom stats quantile
#' @importFrom utils stack
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid unit.pmax
plotPheno <- function(con, comp.sel = FALSE, marker.ord = NA, relevant.only = FALSE, 
  plot.cutoffs = FALSE, delta.panel = TRUE, max.delta.bars = 6, ask = TRUE, grobs = FALSE) {
  
  if (is.null(con$merge)) 
    stop("use contourMerge() first to calculate marker separation")
  
  # comp.sel=TRUE plots data for all outmost contours (root) and outer contours
  # containing inner subsets in a loop
  if (comp.sel && !is.numeric(comp.sel) && !is.null(con$comb)) {
    outer.sq <- unlist(con$comb$sel$outer)
    GGL <- vector("list", length(outer.sq))
    gcount <- 1
    graphics::par(ask = ask)
    for (o in outer.sq) {
      GGL[[gcount]] <- plotPheno(con, comp.sel = o, marker.ord = marker.ord, 
        relevant.only = relevant.only, plot.cutoffs = plot.cutoffs, delta.panel = delta.panel, 
        max.delta.bars = max.delta.bars, ask = ask, grobs = grobs)
      gcount <- gcount + 1
    }
    graphics::par(ask = FALSE)
    if (grobs) 
      return(GGL) else return(invisible(NULL))
  }
  
  title <- "Marker Separation"
  
  # entry for loop or single contour selected
  if (is.numeric(comp.sel) && (length(comp.sel) == 1) && !is.null(con$comb)) {
    if (is.na(comp.sel)) {
      i <- which(is.na(con$comb$sel$outer))
      title <- "Marker Separation (between root contours)"
      
    } else {
      i <- which(con$comb$sel$outer == comp.sel)
      title <- paste0("Marker Separation (within ", rownames(con$comb$sel)[i], 
        ")")
    }
    if (!length(i)) 
      stop("at least 2 contours or an outer contour need to be selected")
    sel <- con$comb$sel$inner[[i]]
  } else if (is.numeric(comp.sel)) {
    sel <- comp.sel
    delta.panel <- FALSE
  } else if (!comp.sel) {
    sel <- which(!grepl("\\*", rownames(con$res)))
    delta.panel <- FALSE
  }
  
  # expression data
  data <- con$comb$data
  
  if (relevant.only) {
    if (length(con$comb$res) > 1) {
      keep <- apply(do.call("rbind", sapply(con$comb$res, "[", 2)), 2, any)
    } else {
      keep <- apply(con$comb$res[[1]]$keep, 2, any)
    }
    
  } else keep <- !logical(ncol(data))
  
  
  if (!is.null(con$beta)) {
    beta <- con$beta
    # message('Notches show range of expression values, \nhinges correspond to the
    # ', round(100 * beta),'th and ', round(100 * (1 - beta)),'th percentiles (beta,
    # 1-beta).')
  } else {
    beta <- 0.25
    # message('Notches show range of expression values, \nhinges default to the 25th
    # and 75th percentiles.')
  }
  
  if (plot.cutoffs && !is.null(con$cutoffs)) {
    beta <- con$beta
    cutoffs <- con$cutoffs  # manual.cutoffs
    cutoffs.chnl <- names(cutoffs)
    n.co <- max(c(2, sapply(cutoffs, length)))  # max number of cutoffs
    
    # if (!all(cutoffs.chnl %in% colnames(data))) { # TODO: provide fs.desc or do
    # outside function?  names(cutoffs) <- switch.names(as.matrix(fs.desc[,1:3]),
    # cutoffs.chnl, colnames(data)) }
    data <- con$comb$data
    keep <- colnames(data) %in% names(cutoffs)
    
  }
  
  # update expression data
  data <- data[, keep, drop = FALSE]
  
  # y limits for marker expression (mlim) fixed for all plots
  mlim <- c(min(apply(data, 2, min)), max(apply(data, 2, max)))
  # option for ordering of markers
  if (is.numeric(marker.ord)) {
    ord <- marker.ord[keep]
  } else ord <- order(colnames(data))
  
  
  beta.box <- function(d) {
    box <- grDevices::boxplot.stats(d)
    return(data.frame(ymin = box$stats[1], ymax = box$stats[5], upper = stats::quantile(d, 
      1 - beta, na.rm = TRUE), lower = stats::quantile(d, beta, na.rm = TRUE), middle = box$stats[3]))
  }
  
  if (!is.null(con$delta)) {
    delta <- con$delta
  } else {
    delta.panel <- FALSE
  }
  
  # boxplot data
  dat <- data.frame(check.names = F, check.rows = F)
  for (cl in sel) {
    subs.id <- rownames(con$res)[cl]
    subset <- c(1:nrow(data)) %in% con$res$bin.idx[[cl]]
    dat.tmp <- data.frame(subset = rep(subs.id, ncol(data)), utils::stack(as.data.frame(data[subset, 
      ])))
    dat <- rbind(dat, dat.tmp)
  }
  # set order of markers
  if ((nlevels(dat$ind) == 1) && (levels(dat$ind) == "data[subset, ]")){
  	dat$ind <- factor(colnames(data))
  } else {
   dat$ind <- factor(dat$ind, levels = colnames(data)[ord])
  }
  # levels(dat$subset); unique(as.numeric(dat$subset)) ggplots legend in separate
  # plot
  gl <- ggplot2::ggplot(dat) + ggplot2::geom_boxplot(ggplot2::aes(x = ind, y = values, fill = subset), lwd = 0.25) + 
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) + ggplot2::theme(legend.title = ggplot2::element_blank(), 
    legend.key = ggplot2::element_blank())

  gl <- gg_legend(gl)  # separate plot
  
  # option for manual cutoffs, list of phenotypes
  g1 <- ggplot2::ggplot(dat) + ggplot2::stat_summary(fun.data = beta.box, geom = "boxplot", aes(x = ind, 
    y = values, fill = subset), position = "dodge", lwd = 0.25)

  if (nlevels(dat$ind) > 1) {
  	# to better discriminate between markers ('dodge')
    g1 <- g1 + ggplot2::geom_vline(xintercept = seq(1.5, nlevels(dat$ind) - 0.5, 1), lwd = 0.25, colour = "lightgray")
  }
  
  if (plot.cutoffs && !is.null(con$cutoffs)) {
    dat.cutoffs <- utils::stack(makeDF(cutoffs[colnames(data)], 1:n.co))  # max number of cutoffs
    dat.cutoffs$ind <- factor(dat.cutoffs$ind, levels = colnames(data)[ord])
    
    g1 <- g1 + ggplot2::geom_point(data = dat.cutoffs, ggplot2::aes(x = ind, y = values), color = "blue", 
      shape = "_", size = 5)
  }
  
  # g1 <- g1 + scale_y_continuous(labels = scaleFUN, limits = mlim )
  
  if (!delta.panel) {
    
    g1 <- g1 + ggplot2::theme(plot.title = ggplot2::element_text(size = 9), legend.position = "none", 
      panel.grid.major.y = ggplot2::element_line(size = 0.5, linetype = 3, colour = "lightgray"), 
      panel.grid.minor.y = ggplot2::element_line(size = 0.25, linetype = 3, colour = "lightgray"), 
      panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white", 
        colour = "black"), axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_text(size = 12), 
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"), axis.text.y = ggplot2::element_text(size = 8)) + 
      ggplot2::ggtitle(title) + ggplot2::ylab("MSI")  # 'asinh MSI'
    
    gA <- ggplot2::ggplotGrob(g1)
    GG <- gridExtra::arrangeGrob(gl, gA, ncol = 1, heights = c(0.5, 3))
    
  } else {
    
    # use absolute delta values
    delta.mat <- abs(con$comb$res[[i]]$delta[, keep, drop = FALSE])
    inner <- con$comb$sel$inner[[i]]
    n.bars <- min(nrow(delta.mat), max.delta.bars)
    # y limits for effective marker difference (elim)
    elim <- c(0, max(c(1.2 * delta, max(delta.mat))))
    # elim <- c( max(-7.5, min(c(-1.2 * delta, min(delta.mat)))), min(10, max(c(1.2 *
    # delta, max(delta.mat)))) )
    
    # barplot data
    dat.eff <- data.frame(check.names = F, check.rows = F)
    for (r in seq_len(nrow(delta.mat))) {
      comp.id <- rownames(delta.mat)[r]
      col.1 <- con$comb$res[[i]]$idx[[r]][1]
      col.2 <- con$comb$res[[i]]$idx[[r]][2]
      dat.eff.tmp <- data.frame(comp.id = rep(comp.id, ncol(delta.mat)), utils::stack(delta.mat[r, 
        ]), col.1 = rep(col.1, ncol(delta.mat)), col.2 = rep(col.2, ncol(delta.mat)))
      dat.eff <- rbind(dat.eff, dat.eff.tmp)
    }
    # set maximum number of delta bars to be displayed per marker
    if (n.bars != nrow(delta.mat)) {
      dat.eff.max <- by(dat.eff, dat.eff$ind, function(x) {
        x[order(abs(x$value), decreasing = TRUE)[1:n.bars], ]
      })
      dat.eff <- do.call("rbind", dat.eff.max)
    }
    # set levels for coloring!
    dat.eff$col.1 <- factor(dat.eff$col.1, levels = (inner))
    dat.eff$col.2 <- factor(dat.eff$col.2, levels = (inner))
    con.cols <- gg_color_hue(length((inner)))
    names(con.cols) <- (inner)
    # set order of markers
    dat.eff$ind <- factor(dat.eff$ind, levels = colnames(data)[ord])
    
    # add theme to boxplot
    g1 <- g1 + ggplot2::theme(plot.title = ggplot2::element_text(size = 9), panel.grid.major.y = ggplot2::element_line(size = 0.5, 
      linetype = 3, colour = "lightgray"), panel.grid.minor.y = ggplot2::element_line(size = 0.25, 
      linetype = 3, colour = "lightgray"), panel.grid.major.x = ggplot2::element_blank(), 
      panel.background = ggplot2::element_rect(fill = "white", colour = "black"), axis.title.y = ggplot2::element_text(size = 12), 
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", 
      axis.title.x = ggplot2::element_blank()) + ggplot2::ggtitle(title) + ggplot2::ylab("MSI")  # 'asinh MSI'
    
    # barplot
    g2 <- ggplot2::ggplot(dat.eff) + ggplot2::geom_abline(slope = 0, intercept = delta, col = "coral", lwd = 0.25) +
    ggplot2::geom_bar(ggplot2::aes(x = ind, y = values, group = comp.id, fill = col.1), position = "dodge", 
      stat = "identity") + # overplot with halved values to get colour code for each comparison
    ggplot2::geom_bar(ggplot2::aes(x = ind, y = values/2, group = comp.id, fill = col.2), position = "dodge", 
      stat = "identity") + ggplot2::scale_fill_manual(values = con.cols)

    if (nlevels(dat$ind) > 1) {
      g2 <- g2 + ggplot2::geom_vline(xintercept = seq(1.5, nlevels(dat$ind) - 0.5, 1), lwd = 0.25, colour = "lightgray")
    }   
     g2 <- g2 + ggplot2::scale_y_continuous(labels = scaleFUN, 
      limits = elim) + ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_line(size = 0.5, 
      linetype = 3, colour = "lightgray"), panel.grid.minor.y = ggplot2::element_line(size = 0.25, 
      linetype = 3, colour = "lightgray"), panel.grid.major.x = ggplot2::element_blank(), 
      panel.background = ggplot2::element_rect(fill = "white", colour = "black"), axis.title.x = ggplot2::element_blank(), 
      axis.title.y = ggplot2::element_text(size = 12), axis.text.x = ggplot2::element_text(angle = 45, 
        hjust = 1, face = "bold"), axis.text.y = ggplot2::element_text(size = 8)) + 
      ggplot2::ylab(expression(Delta))
    
    # https://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
    gA <- ggplot2::ggplotGrob(g1)
    gB <- ggplot2::ggplotGrob(g2)
    maxWidth <- grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
    gA$widths[2:5] <- as.list(maxWidth)
    gB$widths[2:5] <- as.list(maxWidth)
    GG <- gridExtra::arrangeGrob(gl, gA, gB, ncol = 1, heights = c(0.5, 3, 3))
    
  }  # end delta.panel
  
  if (grobs) {
    return(GG)
  } else {
    # plot it
    gridExtra::grid.arrange(GG)
  }
}

## @description helper function to extract legend from a ggplot and plot it separately
## @param g ggplot2 object containing legend
## @return ggplot2 object containing just a legend
## @rdname plotPheno
#' @importFrom ggplot2 ggplot_gtable ggplot_build
gg_legend <- function(g) {
# https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(g))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## @description helper function
## @param list list of numeric vectors
## @param names vector of list names
## @return data.frame
## @rdname plotPheno
makeDF <- function(list, names) {
  m <- (vapply(list, FUN = function(x) unlist(x)[names], FUN.VALUE = numeric(length(names))))
  as.data.frame(m, check.names = FALSE)
}

## @description helper function
## @param x numeric vector
## @return character vector with formatted strings
## @rdname plotPheno
scaleFUN <- function(x) sprintf("%.1f", x)
