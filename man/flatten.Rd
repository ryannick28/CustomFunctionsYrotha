\name{flatten}
\alias{flatten}
\title{Flatten a recursive list}
\description{
This function flattens a recursive list into a non-recursive list, keeping all the elements in the leaves of the original list intact. Therefore, it behaves slightly differently than e.g. unlist. The function was taken from the internet, see: https://stackoverflow.com/questions/47603578/flatten-recursive-list
}
\usage{
flatten(x)
}
\arguments{
\item{x}{An R object, intended to be a recursive list.}
}
\value{
Returns a non-recursive list, listing all elements in the leaves of x.
}
\examples{
X <- list(a = list(b = pi, c = iris[1:3,]), d = 1:10)
flatten(X)
}

