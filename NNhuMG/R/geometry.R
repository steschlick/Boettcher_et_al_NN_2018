#' @title Polygonal geometry
#' @description Internal functions to work with polygons
#' @param outer polygon which is a list with elements x and y containing x- and y-coordinates
#' @param inner list of inner polygon-lists with elements x and y containing x- and y-coordinates
#' @param poly list with elements x and y containing x- and y-coordinates
#' @param points numeric matrix with points in rows and coordinates in first two columns
## @return A polygon is a list with elements x and y containing x- and y-coordinates
#' @details \code{innerPolyExclude} combines polygons so that inner polygons form an exclusion inside an outer polygon.
#'          \code{closest} finds the nearest neighbour of a given point from a set of points in 2D.
#'          The functions \code{arpoly}, \code{inpoly}, and \code{inpolyl} are wrappers around splancs' \code{\link{areapl}},
#'          \code{\link{inpip}}, and \code{\link{inout}}, respectively.
#'          \code{inner} usually is a list of contours, both \code{outer} and \code{poly} a contour, as returned by \code{\link{contourLines}}. 
#' @rdname geometry
#' @keywords internal
innerPolyExclude <- function(outer, inner) {
  outer.pts <- do.call("cbind", outer[[1]])
  # open inner polygon(s)
  inner <- lapply(inner, function(i) list(x = i$x[-1], y = i$y[-1]))
  
  while (length(inner) > 0) {
    dm <- sapply(inner, function(cl) {
      apply(outer.pts, 1, function(x) {
        closest(do.call("cbind", cl), x, TRUE)
      })
    }, simplify = TRUE)
    ai <- arrayInd(which.min(dm), .dim = dim(dm))
    oi <- ai[, 1]  # open outer index
    wi <- ai[, 2]  # which inner to insert
    insert.pts <- do.call("cbind", inner[[wi]])
    i <- closest(insert.pts, outer.pts[oi, ], FALSE)
    ni <- nrow(insert.pts)
    insert.pts <- rbind(insert.pts[i:ni, ], insert.pts[1:i, ])
    # insert inner points
    outer.pts <- rbind(outer.pts[1:oi, ], insert.pts, 
                       outer.pts[oi:nrow(outer.pts), ])
    inner <- inner[-wi]
  }
  return(list(x = outer.pts[, 1], y = outer.pts[, 2]))
}

## @title helper function
## @description Calculates the area of a polygon.
## @param poly list with elements x and y containing x- and y-coordinates
## @return area
## @details This function is a wrapper around splancs' \code{\link{areapl}}.
## It takes the coordinates returned by \code{\link{contourLines}} as input
#' @rdname geometry
#' @importFrom splancs as.points areapl
#' @keywords internal
arpoly <- function(poly) {
  poly <- splancs::as.points(poly$x, poly$y)
  splancs::areapl(poly)
}

## @title helper function
## @description finds points inside a polygon and returns their indices
## @param points numeric matrix with points in rows and coordinates in first two columns
## @param poly list with elements x and y containing x- and y-coordinates
## @return numeric vector
## @details This function is a wrapper around splancs' \code{\link{inpip}}.
#' @rdname geometry
#' @importFrom splancs as.points inpip
#' @keywords internal
inpoly <- function(points, poly) {
  points <- splancs::as.points(matrix(points[, 1:2], ncol = 2))
  poly <- splancs::as.points(poly$x, poly$y)
  splancs::inpip(points, poly, bound = TRUE)
}

## @title helper function
## @description finds points inside a polygon and returns a logical vector
## @param point vector of length 2
## @param poly list with elements x and y containing x- and y-coordinates
## @return logical vector
## @details This function is a wrapper around splancs' \code{\link{inout}}.
#' @seealso 
#'  \code{\link[splancs]{areapl}}
#'  \code{\link[splancs]{inout}}
#'  \code{\link[splancs]{inpip}}
#' @rdname geometry
#' @importFrom splancs as.points inout
#' @keywords internal
inpolyl <- function(points, poly) {
  points <- splancs::as.points(matrix(points[, 1:2], ncol = 2))
  poly <- splancs::as.points(poly$x, poly$y)
  splancs::inout(points, poly, bound = TRUE)
}

## @title helper function
## @description Finds the nearest neighbour of a given point from a set of points in 2D.
#' @param query set of points, a matrix with 2 columns
#' @param point vector of length 2
#' @param distance logical whether to return the shortest distance. If FALSE, the nearest neighbour index is returned, Default: TRUE
## @return distance to or index of the 1st nearest neighbour
#' @rdname geometry
#' @keywords internal
closest <- function(query, point, distance = TRUE) {
  d <- sqrt(rowSums(sweep(query, 2, point)^2))
  i <- which.min(d)
  if (distance) 
    d[i] else i
}
