\name{parReset}
\alias{parReset}
\title{Reset to Old par() Settings}
\description{
Used to reset to previously saved par() settings, ideally stored with parSave()
}
\usage{
parReset(parOld, split = FALSE)
}
\arguments{
\item{parOld}{Previously stored output of par(). Ideally stored with parSave()}

\item{split}{If set to TRUE the plot region is automatically splitted in two rows and two columns}
}
\value{
The output of par() with some values removed in order to make falling back to these settings easier. However, some entries in the output list might still give problems.
}
\examples{
olpar <- parSave()
parReset(olpar)
}

