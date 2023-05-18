\name{flatten}
\alias{flatten}
\title{Flatten a recursive list}
\description{
This function flattens a recursive list into a non-recursive list, keeping all the elements (which do not inherit class "list") in the leaves of the original list intact. Therefore, data frames remain for example unchanged which is a different behavior compared to unlist. The function was taken from the internet, see: https://stackoverflow.com/questions/47603578/flatten-recursive-list
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

