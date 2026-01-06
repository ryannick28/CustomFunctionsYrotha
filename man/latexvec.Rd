\name{latexvec}
\alias{latexvec}
\title{Format a vector for listing its elements in latex text}
\description{
The function formats a vector so that it can be used in a Sweave file e.g. with the \\Sexpr function. The elements of the vector will be separated by specific signs, e.g. commas.
}

\usage{
latexvec(x, colps = ", ", fincolps = " and ", usetxttt = TRUE, fxundscr = TRUE)
}

\arguments{
\item{x}{A vector with at least two elements. It will be turned into a character vector inside the function.}
\item{colps}{A character string which will be used to separate all the elements of x except the last one.}
\item{fincolps}{A character string which will separate the last element of x from the rest.}
\item{usetxttt}{A logical value which determines whether the element names should be written in latex font in the final text document (using the \\texttt function).}
\item{fxundscr}{A logical value which determines whether underscores should be fixed in the elements. Fixing means that escape characters will be inserted before all underscores so that they are correctly printed in the resulting text document.}
}

\value{
Returns a character string which contains the listed elements ready to be interpreted as latex code.
}

\examples{
vc = c("aa", "bb", "cc")
latexvec(x = vc)
}

