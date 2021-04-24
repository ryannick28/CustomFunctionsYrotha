\name{parMarModify}
\alias{parMarModify}
\title{Modify Plot Margins}
\description{
Increases or decreases the margins of the current plot by a specified factor.
}
\usage{
parMarModify(f = 0.3, reduce = FALSE)
}
\arguments{
\item{f}{The factor by which the margins are changed.}
\item{reduce}{A logical value indicating whether the margins should be increased or reduced. Default increases the margins.}
}
\examples{
plot(1:5)
parMarModify()
plot(1:5)
}

