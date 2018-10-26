#' @title Boolean Gate
#' @description Constructs a Boolean \code{flowCore} \code{\link[flowCore]{filter}}
#' @param expr an expression text, see example
#' @param ... further args passed to constructor
#' @param filterId gate name, Default: 'defaultBooleanFilter'
#' @return boolean gate which can be \code{\link[flowWorkspace]{add}ed} to a \code{\link[flowWorkspace]{GatingSet}}
#' @details Constructs an expression for a Boolean gate to be used in a \code{\link[flowWorkspace]{GatingSet}}.
#' @examples 
#' \dontrun{
#' exp <- "/CD3+ & /CD4- & /CD8+" # 3 pops must be present in gs
#' bf <- boolFilter(exp, filterId = "CD8Tcells")
#' # add(gs, bf) # add to GatingSet
#' }
#' @rdname boolFilter
#' @export 
#' @importFrom methods new
boolFilter <- function(expr, ..., filterId = "defaultBooleanFilter") {
  if (missing(filterId)) {
    filterId <- "defaultBooleanFilter"
    }
  methods::new("booleanFilter", filterId = filterId, expr = parse(text = paste0("`", expr, 
    "`")), args = list(...), deparse = expr)
}
