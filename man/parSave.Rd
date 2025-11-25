\name{parSave}
\alias{parSave}
\title{Save par() Settings}
\description{
Returns the current par() settings. Can be used in combination with parReset() to fall back to old par() settings later.
}
\usage{
parSave()
}
\value{
The output of par() with some values removed in order to make falling back to these settings easier. However, some entries in the output list might still give problems.
}
\examples{
olpar <- parSave()
parReset(olpar)
}
