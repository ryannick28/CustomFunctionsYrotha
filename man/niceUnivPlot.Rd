\name{niceUnivPlot}
\alias{niceUnivPlot}
\title{Plot Univariate Variable Nicely}
\description{
Creates a violin plot of a numeric variable and marks the mean value. Can split the plotted variable along a specified factor.
}
\usage{
niceUnivPlot(numVar, catVar = NULL, pairedVar = NULL, violin = TRUE, showMean = TRUE,
bw = 'nrd0', jitFactor = 0.2, add.ylim = 0, ylim.cust = NULL, densScl = 0.5,
main = NULL, sigTest = FALSE, multCmp = FALSE, lCol = NULL)
}
\arguments{
  \item{numVar}{The numeric vector to be plotted}

  \item{catVar}{An optional factor in case numVar should be plotted separately for each factor level}

  \item{pairedVar}{A factor indicating a grouping variable for repeated measures (e.g. person ID). If provided together with a catVar, lines will be plotted indicating the profile of individual cases.}

  \item{violin}{A logical value indicating whether the violin plot should be drawn}

  \item{showMean}{A logical value indicating whether the mean value should be indicated in plot}

  \item{bw}{The smoothing bandwidth to be used. See density() function.}

  \item{jitFactor}{The jittering factor applied.}

  \item{add.ylim}{A number indicating the increase of the y-limits in plot.}

  \item{ylim.cust}{Optionally the y-limits can be set manually. Has to be a numeric vector of length two.}

  \item{densScl}{An optional factor to increase/decrease the width of the violin plots.}

  \item{main}{A character string defining the plot title.}

  \item{sigTest}{A logical value indicating whether group comparisons should be performed (unpaired Wilcoxon tests). The results will be drawn into the plot.}

  \item{multCmp}{A logical value indicating whether the pvalues of the performed group comparisons should be Bonferroni corrected.}

  \item{lCol}{Colour of the lines to be plotted for individual cases.}

}
\examples{
niceUnivPlot(iris$Sepal.Length)
niceUnivPlot(iris$Sepal.Length, catVar = iris$Species)
}
