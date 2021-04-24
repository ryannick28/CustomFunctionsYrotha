\name{asArguments}
\alias{asArguments}
\title{Automatically store the argument-values of a function}
\description{
By copying the arguments of a function and supplying it to asArguments() as a string, each argument is stored with its corresponding value. Important: The function tries to separate arguments by commas followed by a space. Therefore, while commas separating arguments must be followed by a space, all other commas in the string (e.g. in a vector) may not be followed by a space (see examples). Also, the complete string should be flanked by double quotes while strings inside the arguments should be flanked by single quotes to prevent erros. Arguments in the string which are not assigned a value are ignored. The stored arguments can be used to run the expressions inside a function, e.g. for debugging.
}
\usage{
asArguments(argString)
}
\arguments{
\item{argString}{A string containing the arguments which should be stored.}
}
\examples{
argstr <- "a=10, b=c('a','b','c'), c = TRUE"
asArguments(argstr)
list(a, b, c)
}

