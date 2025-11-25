\name{asArguments}
\alias{asArguments}
\title{Store the arguments of a function in the global environment}
\description{
By copying the arguments of a function and supplying them to asArguments, each argument is stored with its corresponding value in the global environment. The stored arguments can be used to run the expressions inside a function, e.g. for debugging.
}
\usage{
asArguments(...)
}
\arguments{
\item{...}{The named arguments which will be store in the global environment.}
}
\examples{
asArguments(a=1, b=1:10, d=list('a', 10, TRUE), gg)
ls()
}

