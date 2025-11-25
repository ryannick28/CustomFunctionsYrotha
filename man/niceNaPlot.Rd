\name{niceNaPlot}
\alias{niceNaPlot}
\title{Plot the amount of missing data in a data set}
\description{
A function to plot the amount of missing data in a data set. The function plots a grid of observations and variables, indicating the missing data by color. To make the plot more easily readable the function orders the observations according to similar NA-patterns using hierarchical clustering. If the data comprises a large amount of unique NA-patterns, the function will first aggregate the data into partitional clusters using kmeans clustering. Afterwards, hierarchical clustering will be applied within and between clusters to determine a good ordering of the observations. The reason for this procedure is to prevent R from crashing when applying hierarchical clustering directly to large data.
}
\usage{
niceNaPlot(x, IDvar = NULL, show_xlab = TRUE, forceUnaggr = FALSE,
nClust = 100, critVal = 10000, verbose = TRUE)
}
\arguments{
\item{x}{A data frame or matrix.}
\item{IDvar}{Optional character string to indicate the name of the observation-ID-variable in x (see examples). The names in IDvar will be used as labels in the generated plot.}
\item{show_xlab}{A logical value indicating whether the xlabels are shown in the generated plot.}
\item{forceUnaggr}{A logical value indicating whether a direct, unaggregated clustering of the NA-patterns should be performed even when there is a number of unique NA-patterns larger than critVal.}
\item{nClust}{How many clusters to find when performing kmeans clustering of the NA-patterns. Only relevant if more unique NA-patterns than critVal.}
\item{critVal}{Number defining how many unique NA-patterns are allowed before switching to an aggregated version of the NA-pattern ordering.}
\item{verbose}{Logical indicating if a progress bar for the kmeans clustering should be shown.}
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

