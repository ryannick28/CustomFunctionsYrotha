\name{niceNaPlot}
\alias{niceNaPlot}
\title{Plot the amount of missing data in a data set}
\description{
A function to plot the amount of missing data in a data set. The function plots a grid of observations and variables, indicating the missing data by color.
}
\usage{
niceNaPlot(x, IDvar = NULL)
}
\arguments{
\item{x}{A data frame or matrix.}
\item{IDvar}{Optional character string to indicate the name of the observation-ID-variable in x (see examples). The names in IDvar will be used as labels in the generated plot.}
}
\value{
If the function call is assigned to a variable, it will store the (ordered) data with the missing values indicated.
}
\examples{
d <-as.matrix(iris)
n1 <- sample(length(d), size=40)
d[n1] <- NA
niceNaPlot(x = d)
### Applied to data frame with IDvar:
d <- as.data.frame(d)
d$ID <- paste0('S', 1:nrow(d))
res <- niceNaPlot(x = d, IDvar = 'ID')
res
}

