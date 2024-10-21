\name{niceUnivPlot}
\alias{niceUnivPlot}
\title{Plot Univariate Variable Nicely}
\description{
Creates a violin plot of a numeric variable and marks the mean value. Can split the plotted variable along a specified factor.
}
\usage{
niceUnivPlot(numVar, catVar = NULL, pairedVar = NULL, violin = TRUE,
fxdCol = NULL, showMean = TRUE, plot.points = TRUE, bw = 'nrd0',
jitFactor = 0.2, add.ylim = 0, ylim.cust = NULL, xlab = NULL,
ylab = NULL, densScl = 0.5, main = NULL, sigGroup = FALSE, sigMu = NULL,
multCmp = FALSE, pairCol = NULL, add.lgnd = FALSE, add = FALSE, lnk.means = NULL,
lnk.means.lwd = 2, pair.lwd = 2, point.ColPalette = NULL, bgalph = 0, ...)
}
\arguments{
  \item{numVar}{The numeric vector to be plotted. numVar can also be supplied as a table (matrix or data.frame). In that case every column is plotted separately.}

  \item{catVar}{An optional factor in case numVar should be plotted separately for each factor level.}

  \item{pairedVar}{A factor indicating a grouping variable for repeated measures (e.g. person ID). If provided together with a catVar, lines will be plotted indicating the profile of individual cases. Additionally, performed group comparisons will switch to paired Wilcoxon tests. If numVar is supplied as a table then pairedVar should be supplied as a factor with same length as the number of rows in numVar (since a wide format is assumed in that case).}

  \item{violin}{A logical value indicating whether the violin plot should be drawn.}

  \item{fxdCol}{One colour in which to plot the points, densities and mean lines.}

  \item{showMean}{A logical value indicating whether the mean value should be indicated in plot.}

  \item{plot.points}{A logical value indicating whether the points should be plotted.}

  \item{bw}{The smoothing bandwidth to be used. See density() function.}

  \item{jitFactor}{The jittering factor applied.}

  \item{add.ylim}{A number indicating the increase of the y-limits in plot.}

  \item{ylim.cust}{Optionally the y-limits can be set manually. Has to be a numeric vector of length two.}

  \item{xlab}{Set the xlabel.}

  \item{ylab}{Set the ylabel.}

  \item{densScl}{An optional factor to increase/decrease the width of the violin plots.}

  \item{main}{A character string defining the plot title.}

  \item{sigGroup}{A logical value indicating whether group comparisons should be performed (two-samples Wilcoxon tests). The results will be drawn into the plot.}

  \item{sigMu}{A numerical value indicating a location mu against which each plotted group will be statistically compared to (one-sample Wilcoxon test). The results will be drawn into the plot.}

  \item{multCmp}{A logical value indicating whether the pvalues of the performed statistical tests should be Bonferroni corrected. The correction is done independently for group comparisons and for mu comparisons.}

  \item{pairCol}{Colour of the lines to be plotted for individual cases (when pairedVar is supplied). As a default every case gets an individual colour (starts repeating at some point). Alternatively, pairCol can be a single colour e.g. supplied as a string. pairCol can also be a factor or numeric with same length as the number of observations. If a factor, then the lines are coloured according to pairCol's levels. If numeric, the lines are coloured following a red-blue gradient according to the values of pairCol (higher values are blue).}

  \item{add.lgnd}{A logical value indicating whether a legend should be added.}

  \item{add}{A logical value indicating whether the content should be plotted into an existing plot.}

  \item{lnk.means}{A colour value. If supplied, the mean values of the groups will be connected by a line in this colour.}

  \item{lnk.means.lwd}{Line width of lines linking the mean values}
  \item{pair.lwd}{Line width of lines linking the individual cases}
  \item{point.ColPalette}{Colour palette that will be used to plot the points (and the violin and mean lines). As a default the hcl.points "dynamic" palette is used with as many colours as there are levels of catVar (or one colour if catVar is not supplied). If point.ColPalette is supplied then it will be repeated by the number of levels of catvar to make sure that there are enough colours for all catvar levels.}
  \item{bgalph}{A transparency factor used for the background of plotted points (only relevant if a pch symbol with a background argument is chosen in the ... argument).}
  \item{...}{Arguments to be passed to the points() call inside the function, where the points are plotted.}
}
\examples{
niceUnivPlot(iris$Sepal.Length)
niceUnivPlot(iris$Sepal.Length, catVar = iris$Species)
d <- iris[1:10,-5]
niceUnivPlot(d, pairedVar = factor(1:nrow(d)), pairCol = 1:nrow(d))
}
