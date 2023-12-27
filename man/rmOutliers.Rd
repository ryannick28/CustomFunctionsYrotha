\name{rmOutliers}
\alias{rmOutliers}
\title{Remove outliers from numeric variables in data frame}
\description{
The function removes outliers from all numeric variables (also integers) in a data frame while leaving all other variables unchanged. Outliers are (columnwise) defined through the criterion used in the boxplot function.
}

\usage{
rmOutliers(x, chngName = TRUE)
}

\arguments{
\item{x}{A data frame.}
\item{chngName}{A logical value indicating whether a tag should be added to the column names of all numeric variables.}
}

\value{
Returns the data frame with no outliers.
}

\examples{
d <- iris[1:10,]
d_no <- rmOutliers(d)
d_no
}

