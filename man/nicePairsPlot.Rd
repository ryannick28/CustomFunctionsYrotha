\name{nicePairsPlot}
\alias{nicePairsPlot}
\title{Plot a Nice Pairs Plot}
\description{
A selfmade version of a pairs plot (not based on pairs()). The upper panel normally shows scatterplots while the lower panel shows significance tests of the pairwise relations. The diagonal panel shows histograms or barplots depending on the variable type. The significance tests are also adpated depending on variable type. Categorical variables are, if not specified differently, grouped together and indicated with shaded backgrounds in the plots.
}
\usage{
nicePairsPlot(x, catVar = NULL, breaks = "Sturges", density = FALSE,
jitter = FALSE, jitFactor = 1, loess = FALSE, swtchPan = FALSE,
exclude = c("none", "numeric", "factor"), keepOrder = FALSE, facsAtBegin = FALSE,
cex.diag = 2, cex.offdiag = 1, cex.mean = 2, cex.lvls = 1)
}
\arguments{
  \item{x}{A dataframe or matrix.}

  \item{catVar}{A factor determining the colour of the points in the scatterplots.}

  \item{breaks}{"breaks" option for hist() function of diagonal panel.}

  \item{density}{A logical value indicating whether density lines should be drawn in the histograms in the diagonal panel.}

  \item{jitter}{A logical value indicating whether a jitter factor should be applied to the scatterplots of numerical variables. For categorical variables there is always some amount of jitter applied.}

  \item{jitFactor}{The jittering factor applied (is slightly decreased for categorical variables).}

  \item{loess}{A logical value indicating whether a loess fit should be plotted in the scatterplots of numerical variables.}

  \item{swtchPan}{A logical value indicating whether the upper and lower panels should be switched.}

  \item{exclude}{Can either be the name(s) of variables to exclude, or a general exclusion option. The general exclusion option can be used to exclude all categorical or numerical variables.}

  \item{keepOrder}{A logical value. Normally, categorical variables are grouped together, this argument can be used to retain the original column order.}

  \item{facsAtBegin}{A logical value. Normally, categorical variables are grouped together at the end of the plot, this argument can be used to move them to the beginning.}

  \item{cex.diag}{Size expansion factor for the titles in the diagonal figures.}

  \item{cex.offdiag}{An expansion factor which will be multiplied with the textsizes of the coefficients displayed in the off-diagonal. Can be used to incerease or decrease the size of the displayed coefficients.}

  \item{cex.mean}{An expansion factor to set the size of the indicated mean values in the scatter plots showing relations between a factor and a numeric variable.}

  \item{cex.lvls}{An expansion factor to set the size of the names of the levels in barplots.}

}
\examples{
nicePairsPlot(iris, density = TRUE, loess = TRUE)
}
