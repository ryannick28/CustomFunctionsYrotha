\name{standizDat}
\alias{standizDat}
\title{Standardize numeric variables in data frame}
\description{
The function standardizes all numerical variables (also integers) in a data frame while leaving all other variables unchanged.
}
\usage{
standizDat(x, chngName = TRUE, onlyCenter = FALSE)
}
\arguments{
\item{x}{A data frame.}
\item{chngName}{A logical value indicating whether a tag should be added to the column names of all standardized variables.}
\item{onlyCenter}{A logical value. If TRUE, the numerical variables are only centered and not scaled.}
}
\value{
Returns the standardized data frame.
}
\examples{
d <- iris[1:10,]
d_st <- standizDat(d)
d_st
}

