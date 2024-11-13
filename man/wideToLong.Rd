\name{wideToLong}
\alias{wideToLong}
\title{Turn data in wide format to long format}
\description{
A function to turn wide format data to long format data. Can also deal with multiple repeatedly measured variables. The identification of repeatedly measured variables is done through the columnnames in the data. The columnnames of repeated-measures variables must start or end with identifier strings to be recognized as such. The identifier strings for the different measurement occasions (e.g. different timepoints) can be defined through the ind or indCust arguments (see examples).
}
\usage{
wideToLong(x, nRep = NULL, ind = 'T*_', indCust = NULL, repColnm = 'repIdentifier',
            ind_atEnd = FALSE, ignore_unbal = FALSE, verbose = FALSE,
            rmv_linkchar = NULL, shw.mssg = TRUE)
}
\arguments{
\item{x}{A data frame containing the data in wide format}
\item{nRep}{The number of repeated measurement occasions (e.g. how many timepoints in longitudinal study)}
\item{ind}{The identifier strings used to recognize the columns in x containing repeatedly measured variables. Must contain one "*" sign which is a placeholder for the numbers 1 to nRep.}
\item{indCust}{A character vector containing the names of custom identifier strings for repeated measurement occasions.}
\item{repColnm}{The name given to the column in the long format data frame indicating the repeated measurement occasions.}
\item{ind_atEnd}{Logical value. Set to TRUE in case the identifier strings are located at the end of the repeatedly measured variable's names.}
\item{ignore_unbal}{Logical value. In case some repeatedly measured variables are not occuring with the same numbers of repetitions, the function will stop and give a message. Set ignore_unbal to TRUE to still continue and fill the missing data with NA's.}
\item{verbose}{Logical value. Set to TRUE to show a percentage bar of progress.}
\item{rmv_linkchar}{A number which will result in removing rmv_linkchar many characters from the repIdentifier column level names. For example if one has timepoint names like "T1_", "T2_", "T3_" setting rmv_linkchar to 1 will remove the underscore from the level names of the final factor in the long format data.}
\item{shw.mssg}{Logical telling if messages should be printed or not. The function will for example turn non-numeric repeated variables into character columns in the long format data and message about this conversion for each such variable.}
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
d <- data.frame(matrix(1:21, nrow = 3, byrow = TRUE))
colnames(d) <- c('S', 'fix', 'startA', 'middleA', 'finishA', 'startB', 'finishB')
wideToLong(x = d, indCust = c('start', 'middle', 'finish'), ignore_unbal = TRUE)
wideToLong(x = d, indCust = c('A', 'B'), ignore_unbal = TRUE, ind_atEnd = TRUE)
### Mix of different classes and NAs:
x <- data.frame(id=factor(paste0('S', 1:3)), fixv='fix', fix2=as.Date(1:3),
x_t1=1:3, x_t2=c(4.1, 5.1, 6.1), x_t3=7:9,
y_t1=factor(LETTERS[1:3]), y_t2=factor(LETTERS[4:6]),
z_t2=c('a', 'b', 'c'),
w_t2=c(as.Date(20), as.Date(88), as.Date(27)), w_t3=11:13)
wideToLong(x = x, nRep = 4, ind = '_t*', ind_atEnd = TRUE, ignore_unbal = TRUE, rmv_linkchar = 1)
}

