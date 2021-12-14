\name{addImgToPlot}
\alias{addImgToPlot}
\title{Add an image to an existing plot}
\description{
This function can be used to add an image stored in a file to an existing R-plot. The position, size and rotation of the image can be specified.
}
\usage{
addImgToPlot(filename, x = NULL, y = NULL, width = NULL, angle = 0, interpolate = TRUE)
}
\arguments{
\item{filename}{A path to an image file.}

\item{x}{Mid x coordinate for image.}

\item{y}{Mid y coordinate for image.}

\item{width}{Width of image (in x coordinate units).}

\item{angle}{Rotation angle of the image. Rotation is around the image's center.}

\item{interpolate}{Passed to graphics::rasterImage. A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing.}
}


