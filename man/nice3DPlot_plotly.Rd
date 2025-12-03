\name{nice3DPlot_plotly}
\alias{nice3DPlot_plotly}
\title{Plot a Nice 3D Plot (using the plotly package)}
\description{
This function creates a nice plot of 3D-data. It is an adaptation of my original nice3DPlot function, now being based on the plotly package whereas the original version was based on the rgl package. Compared to the original version, I removed some features because I rarely used them.
}

\usage{
nice3DPlot_plotly(X = NULL, plotFit = c('no', 'lin', 'int', 'int2', 'int3'),
catVar = factor(1), pointCol = 1, colRamp = NULL, pointSize = NULL,
pointTrans = 255, planeTrans = 100, axesNames = NULL,
axesLeng = NULL, axeLeng_fac = 1.1, add2obj = NULL, rngs = NULL,
plotlyAxes = FALSE, addgrid = TRUE)
}

\arguments{
  \item{X}{The data to be plotted. This should either be a matrix or a dataframe. The function will use the first three columns for plotting, assigning them the axes x, y, z (in this order). Note: In the original version based on rgl I plotted the y variable on the vertical axis. In this plotly-based version, the z-axis is plotted vertically.}

  \item{plotFit}{Add a surface showing the fit of a linear regression (use 'lin' for model without interaction and 'int' for model with interaction). Default is 'no' for no drawn fit. The options 'int2' and 'int3' are only for cases when a factor is provided with catVar. If such a factor is provided it is included in the displayed model fit. The 'lin' option displays a model in which the factor is included as a fixed factor. The 'int' option displays a model with an interaction between the numeric predictors and the factor of catVar included as a fixed factor. The 'int2' option displays a model in which the interaction terms between the numeric predictors and the factor of catVar are included. The 'int3' option displays a model in which all interaction terms are included.}

  \item{catVar}{An optional factor. Providing such a factor will automatically colour the points according to factor levels. Additionally, the factor will be included when plotting the fit of a linear model with the plotFit argument. See the plotFit argument for details.
  }

  \item{pointCol}{Colour of points.}

  \item{colRamp}{A numeric vector of the same length as there are observations in X. It contains the values which will determine the colour of points along a blue-red colour continuum.}

  \item{pointSize}{Size of points.}

  \item{pointTrans}{Transparency of points.}

  \item{planeTrans}{Transparency of plotted surfaces.}

  \item{axesNames}{An optional vector of length 3, giving the axes names of the x, y and z axis (in this order).}

  \item{axesLeng}{Manually set axes length. Will be extended to a vector of length three assigning the lengths of the x, y and z axis (in this order).}

  \item{axeLeng_fac}{A scaling factor which increases the axes lengths.}

  \item{add2obj}{An existing plotly 3d-plot object into which the data should be plotted.}

  \item{rngs}{Either a number or a list of length three in which the ranges of the visual box of the plotly plot can be defined. If supplied as a single number, then the ranges will be set for each axis as c(-rngs, rngs). If supplied as a list, then the three entries in the list contain the ranges (each one being a vector of length 2) of the x, y and z-directions (in this order). In plotly 3d plots, only points are visualized that lie within this visual box.}

  \item{plotlyAxes}{A logical value indicating whether the default axes used in plotly should be visualized. If some added data lies outside the visual box of the current plot, then a warning will be displayed.}

  \item{addgrid}{A logical value indicating whether a grey grid should be visualized spanning the x and y axis.}
}

\examples{
nice3DPlot_plotly(iris, plotFit = 'lin', catVar = iris$Species)
}
