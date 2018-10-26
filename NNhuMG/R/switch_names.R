#' @title Switch Channel and Marker Names
#' @description convenience function
#' @param m matrix with alternative names in columns
#' @param from vector with names to look up in m
#' @param to either one of c("name", "desc", "pretty.desc") or a vector of strings with some target names
#' @return vector of strings
#' @details columns in \code{m} usually are "name" and "desc" from the \code{\link[flowCore]{parameters}} \code{@Data} slot, i.e.
#'       the channel names and description (marker names) of a flowFrame. "pretty.desc" could be marker names prettified to not contain any 
#'       unwanted characters.
#' @rdname switch.names
#' @export 
switch.names <- function(m, from, to) {
  if ((length(to) == 1) && (to %in% c("name", "desc", "pretty.desc"))) {
    cto <- c("name", "desc", "pretty.desc") %in% to
  } else {
    cto <- apply(m, 2, function(x) all(to %in% x))
  }
  cfro <- apply(m, 2, function(x) all(from %in% x))
  if (!(any(cto) && any(cfro))) 
    stop("no matching name found")
  m[match(from, m[, cfro]), cto]
}
