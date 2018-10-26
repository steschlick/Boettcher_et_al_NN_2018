#' @title Aids in the choice of appropriate bandwidths for kernel density estimation
#' @description Computes kernel density estimates (kde's) for a range of bandwidths, 
#' detects regions of continuously appearing peaks and suggests bandwidths accordingly
#' @param data matrix with tSNE coordinates
#' @param nH number of different bandwidths for which to precompute kde's, Default: 100
#' @param binned logical flag for binned kernel estimation, passed to \code{\link{Hpi}} and \code{\link{kde}}, Default: TRUE
#' @param bgridsize vector of binning grid sizes, passed to \code{\link{Hpi}} and \code{\link{kde}}, Default: c(151, 151)
#' @param store.fhats logical whether to return list of density estimates, Default: TRUE
#' @param MAX.fine maximum number of peaks allowed for "fine" bandwidth suggestion, Default: 35
#' @param MAX.coar maximum number of peaks allowed for "coarse" bandwidth suggestion, Default: 15
#' @param tol.spur number of 'spurious' peaks to be tolerated in "coarse" bandwidth suggestion (see details), Default: 0
#' @param plotit logical whether to plot an overview showing numbers of counted peaks per bandwidth, 
#' and peak coordinates for three suggested bandwidths, Default: TRUE
#' @param ... further args passed to \code{\link{Hpi}}
#' @return A list containing elements: 
#' \describe{
#'   \item{\code{peaks}}{list of peak coordinates}
#'   \item{\code{fhats}}{list of \code{\link{kde}} density estimates}
#'   \item{\code{H}}{plug-in bandwidth matrix as returned by \code{\link{Hpi}}}
#'   \item{\code{Hfactor}}{vector of correction factors}
#'   \item{\code{Hindex}}{vector with indices for \code{H1} (fine), \code{H2} (coarse), and \code{H} (plug-in) bandwidths}
#' }
#' @details 
#'   Although in most statistical inference the bandwidth parameter is determined by the data, 
#'   kernel density estimation using plug-in bandwidth selectors does not guarantee to match 
#'   human-perceived multimodal distributions. For example in very fine-grained tSNE maps,   
#'   small but potentially relevant subsets might be "smoothed away" in the presence of few 
#'   dominant modes. In this context, the bandwidth remains a tunable parameter and the choice
#'   of an appropriate density estimate depends eventually on a "visual criterion". This function 
#'   tries to aid in the decision by performing following steps: 
#' \enumerate{
#'   \item estimate bivariate plug-in bandwidth matrix \code{H} using \code{\link{Hpi}}
#'   \item detect (\code{\link{PeakFind}}) and count number of local maxima over range of 
#'         bandwidth matrices (in multiples of \code{H})
#'   \item remove noisy peaks (via \code{\link{epsKNNee}} and \code{\link{dbscan}}) to improve 
#'         significant feature detection in the next step 
#'   \item obtain significant high negative curvature regions on all detected peaks 
#'         using \code{\link{featureSignif}}
#'   \item count, label and visualize peaks falling into these regions: whereas the 
#'         overall number of peaks (in red) is a monotonically decreasing function 
#'         of the bandwidth, the numbers of continuously appearing local maxima 
#'         (in green) are rather stable. 
#'   \item suggest \code{H1} as the largest bandwidth for which the number of continuous 
#'         peaks is at maximum (usually accompanied by an overall larger number 
#'         of peaks) and \code{H2} as the largest bandwidth for which the number (nearly) converges 
#'         with the overall number of local maxima
#'   \item optionally inspect computed density estimates via \code{\link{shinyBW}} and 
#'         make a judicious choice.
#' }
#' @examples 
#' \dontrun{
#' set.seed(123)
#' x <- rbind(matrix(rnorm(500, sd = 0.3), ncol = 2, 
#'                   dimnames=list(NULL, c("A", "B"))),
#'            matrix(rnorm(36, mean = 1, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = -1, sd = 0.3), ncol = 2))
#' BW <- chooseBandWidth(x, plotit = TRUE)
#' BW[3:5]
#' 
#' # interactive selection of a bandwidth, requires global var 'BW' with slot \code{fhats} present
#' idx <- shinyBW()
#' fhat <- BW$fhats[[idx]]
#' }
#' @seealso 
#'  \code{\link[ks]{Hpi}}
#'  \code{\link[ks]{kde}}
#'  \code{\link[dbscan]{dbscan}}
#'  \code{\link{epsKNNee}}
#'  \code{\link[feature]{featureSignif}}
#'  \code{\link[grDevices]{contourLines}}
#' @references 
#'   K. Shekhar, P. Brodin, M.M. Davis, A.K. Chakraborty (2014), Automatic classification of 
#'   cellular expression by nonlinear stochastic embedding (ACCENSE), 
#'   Proceedings of the National Academy of Sciences, \url{http://www.cellaccense.com/oldver.html},
#'   
#'   T. Duong, A. Cowling, I. Koch, M.P. Wand (2008), Feature significance for multivariate 
#'   kernel density estimation. Computational Statistics and Data Analysis, 52, 4225-4242. 
#' @rdname chooseBandWidth
#' @export 
#' @importFrom ks Hpi kde
#' @importFrom dbscan dbscan
#' @importFrom feature featureSignif
#' @importFrom grDevices contourLines
#' @import graphics
#' @importFrom stats setNames
chooseBandWidth <- function(data, nH = 100, binned=TRUE, bgridsize=c(151, 151), store.fhats = TRUE, 
  MAX.fine = 35, MAX.coar = 15, tol.spur = 0, plotit = TRUE, ...) {
  
  if (nH < 100) {
  	message("provided nH too small, precomputing 100 kde's")
  	nH <- 100
  }
  # use 5 as max number of successive p.n
  K.succ <- 5L
  
  kernel.range <- min(diff(range(data[, 1])), diff(range(data[, 2])))
  
  H <- ks::Hpi(data, binned = binned, bgridsize = bgridsize, ...)
  
  kernel.min <- min(min(diag(H))/10, kernel.range/100)
  kernel.max <- max(max(diag(H)), kernel.range/10)
  bw.range <- seq(kernel.min, kernel.max, length = nH)
  
  # include Hpi and indicate in plot later on
  idx.Hpi <- which.min(abs(bw.range - max(diag(H))))
  
  # test a sequence of bivariate bandwidths
  bw.range.fac <- max(abs(H))/bw.range
  bw.range.fac[idx.Hpi] <- 1
  
  p.l <- f.l <- vector("list", length = nH)
  
  for (i in 1:nH) {
    cat("Counting peaks for bandwidth No.: ", i, "\n")
    fhat <- ks::kde(data, binned = TRUE, bgridsize = bgridsize, H = H/bw.range.fac[i])
    p.l[[i]] <- cbind(PeakFind(fhat), i)
    if (store.fhats) {
      f.l[[i]] <- fhat
    }
  } 
  
  ### routine to find a fine and a coarse bw 
  # keep plug-in bw Hpi for plotting and return
  p.Hpi <- p.l[[idx.Hpi]]
  f.Hpi <- f.l[[idx.Hpi]]
  
  # number of peaks per bandwidths
  p.b <- unlist(sapply(p.l, nrow))
  
  # number of peaks that changes the least over sequence of bw's
  p.modal <- as.numeric(names(which.max(table(p.b))))
  p.modal.ids <- which(p.b == p.modal)
  # truncate at a max number of successive peaks
  p.modal.idx <- p.modal.ids[min(length(p.modal.ids), K.succ)]
  p.l <- p.l[1:p.modal.idx]
  f.l <- f.l[1:p.modal.idx]
  bw.range.fac <- bw.range.fac[1:p.modal.idx]
  
  # number of peaks per bandwidths, truncated
  p.n <- unlist(sapply(p.l, nrow))
  p.r <- do.call("rbind", p.l)[, c(1, 2)]  # peak locations
  p.i <- do.call("rbind", p.l)[, c(4)]  # peak indices
  
  # remove noisy peaks using dbscan outliers
  eps <- epsKNNee(p.r, k = K.succ, FALSE)
  p.dbsc <- dbscan::dbscan(p.r, eps = eps, minPts = K.succ)$cluster
  p.zero <- p.dbsc == 0
  
  # determine significant curvature, i.e. regions of peaks over range of bandwidths
  # without noisy peaks
  p.feat <- feature::featureSignif(p.r[!p.zero, ], signifLevel = 0.05)
  
  # extract polygons, i.e. countour lines
  CL <- grDevices::contourLines(p.feat$fhat[[1]][[1]], p.feat$fhat[[1]][[2]], p.feat$curv, 
    nlevels = 1, levels = 0.5)
  
  if (!length(CL)) {
  	stop("no regions of continuous peaks found. try again with increased nH")
  }
  
  CL.ids <- sapply(CL, function(cl) {
    list(inpoly(p.r, cl))
  })
  
  # assign dbscan p.c to significant curvature ids, include also p.zero if inside a
  # region
  p.sigcurv <- integer(length(p.dbsc))
  # get indices for each dbscan.cluster
  for (i in 1:max(p.dbsc)) {
    idx <- which(sapply(CL.ids, function(cl) any(which(p.dbsc == i) %in% cl)))
    if (!length(idx)) {
      p.sigcurv[p.dbsc == i] <- 0
    } else if (length(idx) == 1) {
      p.sigcurv[p.dbsc == i] <- idx
    } else {
      # assign to closest sigreg if a dbscan.cluster has multiple sigregs
      mat.d <- sapply(CL.ids[idx], function(cl) {
        apply(p.r[p.dbsc == i, ], 1, function(x) closest(p.r[cl, ], x, TRUE))
      })
      p.sigcurv[p.dbsc == i] <- idx[apply(mat.d, 1, which.min)]
    }
  }
  # assign indices for dbsc.zeros separately
  if (any(p.zero)) {
    p.sigcurv[p.zero] <- sapply(which(p.zero), function(z) {
      idx <- which(sapply(CL.ids, function(cl) z %in% cl))
      if (length(idx)) 
        idx else 0
    })
  # count peaks that belong to sigreg and are no noise
    target <- sort(unique(p.sigcurv))[-1]  # w/o 0, i.e. noise
  } else {
    target <- sort(unique(p.sigcurv))	
  }

  p.d <- sapply(min(p.i):max(p.i), function(i) sum(target %in% p.sigcurv[p.i == 
    i]))
  
  # coarse bw index
  p.d.opt.idx.2 <- which(((p.n - p.d) <= tol.spur) & (p.n <= MAX.coar))[1]
  
  # fine bw index
  p.d.opt.ids.1 <- which((p.d == max(p.d[p.n <= MAX.fine])) & (p.n <= MAX.fine))
  p.d.opt.idx.1 <- p.d.opt.ids.1[length(p.d.opt.ids.1)]
  
  if (plotit) {
  	mfrow <- graphics::par()$mfrow
    graphics::layout(matrix(c(1, 1, 1, 2, 3, 4), 2, 3, byrow = TRUE))
    graphics::plot(p.b, col = 8, xlab = "bandwidth", ylab = "# peaks", ylim = range(p.b), 
      xlim = c(1, nH), xaxt = "n")
    at <- floor(seq.int(1, nH, length.out = 6))
    graphics::axis(1, at = at, labels = round(bw.range[at], digits = 1))
    # difference between red and green is number of 'spurious' peaks
    graphics::points(p.d, ylim = range(p.b), col = 3)
    graphics::points(p.n, col = 2)
    graphics::points(c(p.d.opt.idx.2, p.d.opt.idx.1, idx.Hpi), p.b[c(p.d.opt.idx.2, p.d.opt.idx.1, 
      idx.Hpi)], pch = 19, col = 4)
    graphics::text(c(p.d.opt.idx.2, p.d.opt.idx.1, idx.Hpi), p.b[c(p.d.opt.idx.2, p.d.opt.idx.1, 
      idx.Hpi)], labels = c("H2", "H1", "Hpi"), cex = 0.75, adj = c(0.3, -0.6))
    xlim <- range(data[, 1])
    ylim <- range(data[, 2])
    col.sigcurv <- p.sigcurv
    col.sigcurv[col.sigcurv > 0] <- 3
    col.sigcurv[col.sigcurv == 0] <- 2
    graphics::plot(p.r, xlim = xlim, ylim = ylim, col = col.sigcurv, main = "Hpi")
    graphics::contour(p.feat$fhat[[1]][[1]], p.feat$fhat[[1]][[2]], p.feat$curv, nlevels = 1, 
      levels = 0.5, xlim = xlim, ylim = ylim, add = TRUE, drawlabels = FALSE, 
      col = 3)
    graphics::points(p.Hpi[, 1:2], pch = 19, col = 4, cex = 1.5)
    graphics::plot(data, xlim = xlim, ylim = ylim, pch = "*", col = 8, main = "H1: fine")
    graphics::points(p.l[[p.d.opt.idx.1]][, 1:2], pch = 19, col = 4)
    graphics::plot(data, xlim = xlim, ylim = ylim, pch = "*", col = 8, main = "H2: coarse")
    graphics::points(p.l[[p.d.opt.idx.2]][, 1:2], pch = 19, col = 4)
    graphics::points(p.l[[p.d.opt.idx.2]][, 1:2], pch = 19, col = 4)
    graphics::par(mfrow=mfrow)
  }
  
  if (length(p.l) < idx.Hpi) {
    f.l <- c(f.l, list(f.Hpi))
    p.l <- c(p.l, list(p.Hpi))
    bw.range.fac <- c(bw.range.fac, 1)
    idx.Hpi <- length(p.l)
  }
  
  opt.idx <- c(p.d.opt.idx.1, p.d.opt.idx.2, idx.Hpi)
  names(opt.idx) <- c("H1", "H2", "Hpi")
  
  if (!store.fhats) {
    f.l <- NULL
  }
  return(stats::setNames(list(p.l, f.l, H, bw.range.fac, opt.idx), c("peaks", "fhats", 
    "H", "Hfactor", "Hindex")))
}

#' @rdname chooseBandWidth
#' @export 
#' @import rgl
#' @import shiny
shinyBW <- function() {

  ui <- shiny::pageWithSidebar(
    shiny::headerPanel("Choose a bandwidth!"),
    # sidebar with slider input for index of kde
    shiny::sidebarPanel(
      shiny::sliderInput("index", 
                "Bandwidth Index:", 
                min = 1, 
                max = 100, 
                value = 50),
      shiny::actionButton("submit", "Use Selected Bandwidth")
    ),
    # kde 3d plot
    shiny::mainPanel(
         rgl::rglwidgetOutput("plot",  width = 800, height = 600)  
        )
  )

  server <- (function(input, output, session) {
    if (!exists("BW"))
      stop(paste("'BW' doesn't exist. This Shiny App is intended to be run locally",
                   "as a part of an R script."))
    if (is.null(BW$fhats))
       stop("kde's missing. run chooseBandWidth() with 'store.fhats=TRUE'")
    shiny::updateSliderInput(session, "index", value = unname(BW$Hindex[3]),
                             min = 1, max = length(BW$Hfactor), step = 1)
    fhat <- shiny::reactive({ BW$fhats[[input$index]] })
    output$plot <- rgl::renderRglwidget({
      fh <- fhat()
      rgl::rgl.open(useNULL=T)
      rgl::bg3d("white")
      rgl::persp3d(fh$eval.points[[1]], fh$eval.points[[2]], fh$estimate, 
                   xlab="tSNE1", ylab="tSNE2", zlab="", aspect = c(1, 1, .25),
                   alpha=1, col=colorDens(fh$estimate), box = TRUE, axes = TRUE) 
      rgl::rglwidget()
    })
    shiny::observe({
      if (input$submit == 0)
        return()
        shiny::stopApp(input$index)
    })
  })   

 shiny::runApp(list(ui = ui, server = server))

}







