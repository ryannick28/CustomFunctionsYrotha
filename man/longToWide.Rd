\name{longToWide}
\alias{longToWide}
\title{Turn data in long format to wide format}
\description{
A function to turn long format data to wide format data. Can deal with multiple repeatedly measured variables and fixed variables of different types. An observation-ID variable (e.g. subject names), a repeated-measures indicator variable (e.g. "Timepoint") and the names of the repeatedly measured variables have to be provided. The function will give specific warnings in case there are peculiar patterns in the data (e.g. different values of a presumably fixed variable in an ID).
}
\usage{
longToWide(x, IDvar, repColnm, repVars)
}
\arguments{
\item{x}{A data frame containing the data in long format}
\item{IDvar}{A string with the name of the observation ID variable (e.g. subject names)}
\item{repColnm}{The name of the column in the long format data frame indicating the repeated measurement occasions.}
\item{repVars}{A vector including the variable names of the variables with repeated measurements}
}
\value{
The transformed data frame in wide format.
}
\examples{
### Generate artificial data:
d <- data.frame(matrix(1:24, nrow = 3, byrow = TRUE))
colnames(d) <- c('fix', 'S','T1_A', 'T2_A', 'T3_A', 'T1_B', 'T2_B', 'T3_B')
d$S <- paste0('S', 1:3)
dat <- wideToLong(x = d, nRep = 3, repColnm = 'timepoint')
dat$S <- as.factor(dat$S)
dat$B <- factor(dat$B, labels = c('a', 'b','c', 'd', 'e', 'f', 'g', 'h', 'j'))
### Apply function:
longToWide(x = dat, IDvar = 'S', repColnm = 'timepoint', repVars = c('B','A'))
}

