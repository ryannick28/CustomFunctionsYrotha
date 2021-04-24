\name{setUpGraph}
\alias{setUpGraph}
\title{Set up a Window for Plotting}
\description{
Opens a new external window for plotting. Mainly consists of calling dev.new() function. As an invisible output the current par() settings are returned by the use of parSave().
}
\usage{
setUpGraph(row = 2, col = 2)
}
\arguments{
\item{row}{Set number of rows in plot window.}
\item{col}{Set number of columns in plot window.}
}
\value{
Returns (invisibly) the current par() settings by the use of parSave().
}
\examples{
olpar <- setUpGraph(row = 2, col = 2)
plot(1:5)
parReset(olpar)
plot(1:5)
}

