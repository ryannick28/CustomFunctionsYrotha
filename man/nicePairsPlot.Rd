\name{nicePairsPlot}
\alias{nicePairsPlot}
\title{Plot a Nice Pairs Plot}
\description{
Basically a wrapper function for pairs(). The upper panel shows the scatterplots with optional loess fit. The diagonal panel shows a histogramm with optional density line. The lower panel shows the pairwise correlation coefficients with adjusted digit sizes and significant correlations marked with stars.
}
\usage{
nicePairsPlot(x, catVar = NULL, breaks = "Sturges", density = FALSE, jitter = FALSE,
jitFactor = 1, loess = FALSE, swtchPan = FALSE, txtInc=1)
}
\arguments{
  \item{x}{A dataframe or matrix.}

  \item{catVar}{A factor dictating the colour of the points in the scatterplots.}

  \item{breaks}{"breaks" option for hist() function of diagonal panel.}

  \item{density}{A logical value indicating whether density lines should be drawn in the diagonal panel.}

  \item{jitter}{A logical value indicating whether a jitter factor should be applied to the scatterplots. Useful for visualizing distribution of factors.}

  \item{jitFactor}{The jittering factor applied.}

  \item{loess}{A logical value indicating whether a loess fit should be plotted in the scatterplots.}

  \item{swtchPan}{A logical value indicating whether the upper and lower panels should be switched.}

  \item{txtInc}{A number which will be multiplied with the textsizes of the correlation coefficients displayed in the plot. Can be used to incerease or decrease the size of the displayed coefficients.}
}
\examples{
nicePairsPlot(iris, density = TRUE, loess = TRUE)
}
