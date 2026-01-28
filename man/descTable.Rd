\name{descTable}
\alias{descTable}
\title{Create a table of descriptive statistics}
\description{
The function creates a table of descriptive statistics and allows certain modifications to facilitate the use of the table in a latex document (intended to be used with xtable).
}

\usage{
descTable(x, numDesc1 = c("median", "mean"), numDesc2 = c("quartile13", "sd"),
          showMiss = TRUE, boldVarnmsLtx = FALSE, addIndetLtx = FALSE,
          escPercLtx = FALSE)
}

\arguments{
\item{x}{A data table. It will be turned into a data frame within the function using as.data.frame.}
\item{numDesc1}{The main descriptive statistic to be computed for numeric variables.}
\item{numDesc2}{The secondary descriptive statistic to be computed for numeric variables (presented in brackets).}
\item{showMiss}{Logical value indicating whether the numbers of missing observations should be presented.}
\item{boldVarnmsLtx}{A logical value which determines whether the variable names should be printed bold in a latex table.}
\item{addIndetLtx}{A logical value which determines whether there should be some indentation for the level names of factors in a latex table.}
\item{escPercLtx}{A logical value which determines whether for printing in latex the percent signs should be escaped using backward slashes.}
}


\details{
For non-numeric variables, the function automatically turns them into a
factor (using as.factor) and returns the number of cases and relative proportion
for each factor level. For numeric variables, the "quartile13" option returns
the first and third quartiles in bracktes, while the the "sd" option returns
the estimated standard deviation.
}

\value{
Returns a data frame containing the descriptive statistics.
}

\examples{
descTable(x = iris)
### For latex:
descTable(x = iris, boldVarnmsLtx = TRUE, addIndetLtx = TRUE,
          escPercLtx = TRUE)
}

