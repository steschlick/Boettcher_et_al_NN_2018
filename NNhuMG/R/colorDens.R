#' @title A color ramp
#' @description Generates a sequence of MATLAB's jet-like colors
#' @param w numeric vector
#' @return a vector of hex colors
#' @details This function maps the input values to [1, 100] and a applies the colorRamp: c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
#'     "red", "#7F0000").
#' @examples 
#' \dontrun{
#' x <- c(23:234)
#' image(as.matrix(x), col=colorDens(x))
#' }
#' @seealso 
#'  \code{\link[grDevices]{colorRamp}}
#' @rdname colorDens
#' @export 
#' @importFrom grDevices colorRampPalette
colorDens <- function(w) {
  clrs <- grDevices::colorRampPalette(c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
    "red", "#7F0000"))
  v <- w
  if (any(w < 0)) {
    v <- v - min(w[w < 0])
  }
  m <- max(v)
  v <- floor(v * 99/m) + 1
  clrs(100)[v]
}
