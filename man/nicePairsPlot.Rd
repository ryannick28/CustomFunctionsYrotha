\name{nicePairsPlot}
\alias{nicePairsPlot}
\title{Plot a Nice Pairs Plot}
\description{
Basically a wrapper function for pairs(). The lower panel shows the scatterplots with optional loess fit. The diagonal panel shows a histogramm with optional density line. The upper panel shows the pairwise correlation coefficients with adjusted digit sizes and significant correlations marked with stars.
}
\usage{
nicePairsPlot(x, catVar = NULL, breaks = "Sturges", density = FALSE, jitter = FALSE,
jitFactor = 1, loess = FALSE, swtchPan = FALSE)
}
\arguments{
  \item{x}{A dataframe or matrix.}

  \item{catVar}{A factor dictating the colour of the points in the scatterplots.}

  \item{breaks}{"breaks" option for hist() function of diagonal panel.}

  \item{density}{A logical value indicating whether the density lines should be drawn in the diagonal panel.}

  \item{jitter}{A logical value indicating whether a jitter factor should be applied to the scatterplots. Useful for visualizing distribution of factors.}

  \item{jitFactor}{The jittering factor applied.}

  \item{loess}{A logical value indicating whether a loess fit should be plotted in the lower panel.}

  \item{swtchPan}{A logical value indicating whether the upper and lower panels should be switched.}
}
\examples{
nicePairsPlot(iris, density = TRUE, loess = TRUE)
}
