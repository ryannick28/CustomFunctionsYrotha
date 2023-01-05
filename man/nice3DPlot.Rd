\name{nice3DPlot}
\alias{nice3DPlot}
\title{Plot a Nice 3D Plot}
\description{
This function creates a nice plot of 3D-data, including a visualization of the density. If no data is supplied, then only the coordinate system is drawn, which can be used for further plotting. If data is provided, it must be numeric. One can select if only the points, only the density or both shall be plotted. One can also use it to plot simulated data (normally distributed). In that case one has to provide a covariance matrix and a vector of the means.
}

\usage{
nice3DPlot(X = NULL, whatToPlot = c("P","D", "PD"),
plotFit = c('no','lin','int','int2','int3'), catVar = factor(1),
covMat = NULL, means = NULL, nSim = 500, pointCol = 1,
pointSize = NULL, spheres = TRUE, pointTrans = 0.8, Dtransp_fac = 0.06,
colres = 50, gridRes = 30, h = 0.3, DpointSize = 30, axesNames = NULL,
axesTicks = FALSE, gridCol = "grey", axesLeng = NULL, zoom = 1,
add = FALSE, htmlFilename = NULL, ...)
}

\arguments{
  \item{X}{The data to be plotted. This should either be a matrix or a dataframe. The function will use the first three columns for plotting, assigning them the axes x, y, z (in this order).}

  \item{whatToPlot}{A string defining what to plot. One can either plot only the data points ('P'), only the density ('D') or both ('PD').}

  \item{plotFit}{Add a surface showing the fit of a linear regression (use 'lin' for model without interaction and 'int' for model with interaction). Default is 'no' for no drawn fit. The options 'int2' and 'int3' are only for cases when a factor is provided with catVar. If such a factor is provided it is included in the displayed model fit. The 'lin' option displays a model in which the factor is included as a fixed factor. The 'int' option displays a model with an interaction between the numeric predictors and the factor of catVar included as a fixed factor. The 'int2' option displays a model in which the interaction terms between the numeric predictors and the factor of catVar are included. The 'int3' option displays a model in which all interaction terms are included.}

  \item{catVar}{An optional factor. Providing such a factor will automatically colour the points according to factor levels. Additionally, the factor will be included when plotting the fit of a linear model with the plotFit argument. See the plotFit argument for details.
  }

  \item{covMat}{Optional covariance matrix (has to have 3rows/columns, order: x,y,z) to plot simulated data from a multivariate normal distribution.}

  \item{means}{A vector of length three with the x, y and z mean values of the simulated data.}

  \item{nSim}{Number of simulated data points.}

  \item{pointCol}{Colour of points.}

  \item{pointSize}{Size of points.}

  \item{spheres}{Logical indicating whether rgl.spheres() should be used for drawing points (instead of rgl.points()). Depending on hardware (apparently) rgl.points() may display datapoints as boxes. This is why per default spheres is set to TRUE.}

  \item{pointTrans}{Transparency of points.}

  \item{Dtransp_fac}{Determines the transparency of the density pixels. With smaller values the displayed density will be more transparent.}

  \item{colres}{Determines the number of steps in the colour progression of the density pixels.}

  \item{gridRes}{Resolution of the grid for which the density values are calculated.}

  \item{h}{Bandwidth of the 3D-Kernel Density Estimator. Changing the value of h can often lead to drastically improved visualization of the density.}

  \item{DpointSize}{Size of the points of the density grid.}

  \item{axesNames}{An optional vector of length 3, giving the axes names of the x, y and z axis (in this order).}

  \item{axesTicks}{Optional logical value to add ticks to the axes.}

  \item{gridCol}{Color of the floor grid in the 3D-coordinate system.}

  \item{axesLeng}{Manually set axes length. Will be extended to a vector of length three assigning the lengths of the x, y and z axis (in this order).}

  \item{zoom}{How much to zoom into the graph (so far only possible with origin as center). Note: On my machine the field of view can be adapted by moving cursor while having the right-click pressed.}

  \item{add}{Set to TRUE in order to add points to an existing 3D Plot.}

  \item{htmlFilename}{A string giving the name/path of a html file to which the 3D plot will be saved to. If NULL than no html file is created. Dont forget to add the ".html" to the name.}

  \item{...}{Additional arguments which will be passed to the rglwidget() function, which is responsible for creating a html version of the 3D plot.}
}
\details{
The 3D-density is created by generating a grid of points around the data and colouring each point of the grid according to the density value at that position. In addition, the transparency of the gridpoints is set according to the density values.
}
\examples{
# Only coordinate-system:
nice3DPlot()

# Plot data with density:
nice3DPlot(X = data.frame(x=rnorm(100), y=rnorm(100), z=rnorm(100)),
           whatToPlot = 'PD', h = 0.7, DpointSize = 30, zoom=0.9)

# Simulated data:
m <- rep(0,3)
covM <- matrix(c(1,0.5,0.5,0.5, 1,0.5,0.5,0.5,5), ncol=3)
nice3DPlot(covMat = covM, means = m, whatToPlot = 'PD', h = 0.7,
           DpointSize = 30, zoom=1.5, axesLeng = 5, pointTrans = 0.5)
}
