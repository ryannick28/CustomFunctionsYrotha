\name{wideToLong}
\alias{wideToLong}
\title{Turn data in wide format to long format}
\description{
A function to turn wide format data to long format data. Can also deal with multiple repeatedly measured variables. The identification of repeatedly measured variables is done through the columnnames in the data. The columnnames of repeated-measures variables must start with identifier strings to be recognized as such. The identifier strings for the different measurement occasions (e.g. different timepoints) can be defined through the ind or indCust arguments (see examples).
}
\usage{
wideToLong(x, nRep = NULL, ind = 'T*_', indCust = NULL, repColnm = 'repIdentifier')
}
\arguments{
\item{x}{A data frame containing the data in wide format}
\item{nRep}{The number of repeated measurement occasions (e.g. how many timepoints in longitudinal study)}
\item{ind}{The identifier strings used to recognize the columns in x containing repeatedly measured variables. Must contain one "*" sign which is a placeholder for the numbers 1 to nRep.}
\item{indCust}{A character vector containing the names of custom identifier strings for repeated measurement occasions.}
\item{repColnm}{The name given to the column in the long format data frame indicating the repeated measurement occasions.}
}
\value{
The transformed data frame in long format.
}
\examples{
d <- data.frame(matrix(1:24, nrow = 3, byrow = TRUE))
colnames(d) <- c('S', 'fix', 'T1_A', 'T2_A', 'T3_A', 'T1_B', 'T2_B', 'T3_B')
d$S <- paste0('S', 1:3)
wideToLong(x = d, nRep = 3, repColnm = 'timepoint')
### Different indicator:
colnames(d) <- c('S', 'fix', 'occ1A', 'occ2A', 'occ3A', 'occ1B', 'occ2B', 'occ3B')
wideToLong(x = d, nRep = 3, ind = 'occ*')
### Custom indicators:
colnames(d) <- c('S', 'fix', 'startA', 'middleA', 'finishA', 'startB', 'middleB', 'finishB')
wideToLong(x = d, indCust = c('start', 'middle', 'finish'))
}

