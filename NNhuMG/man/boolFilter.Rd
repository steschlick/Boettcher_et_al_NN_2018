% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boolFilter.R
\name{boolFilter}
\alias{boolFilter}
\title{Boolean Gate}
\usage{
boolFilter(expr, ..., filterId = "defaultBooleanFilter")
}
\arguments{
\item{expr}{an expression text, see example}

\item{...}{further args passed to constructor}

\item{filterId}{gate name, Default: 'defaultBooleanFilter'}
}
\value{
boolean gate which can be \code{\link[flowWorkspace]{add}ed} to a \code{\link[flowWorkspace]{GatingSet}}
}
\description{
Constructs a Boolean \code{flowCore} \code{\link[flowCore]{filter}}
}
\details{
Constructs an expression for a Boolean gate to be used in a \code{\link[flowWorkspace]{GatingSet}}.
}
\examples{
\dontrun{
exp <- "/CD3+ & /CD4- & /CD8+" # 3 pops must be present in gs
bf <- boolFilter(exp, filterId = "CD8Tcells")
# add(gs, bf) # add to GatingSet
}
}
