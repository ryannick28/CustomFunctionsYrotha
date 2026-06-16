\name{descTable}
\alias{descTable}
\title{Create a table of descriptive statistics}
\description{
The function creates a table of descriptive statistics and allows certain modifications to facilitate the use of the table in a latex document (intended to be used with xtable).
}

\usage{
descTable(x, indx = list("All observations" = 1:nrow(x)),
                      numDesc1 = c("median", "mean"),
                      numDesc2 = c("quartile13", "sd"), showMiss = TRUE,
                      showZeroMiss = FALSE, boldVarnmsLtx = FALSE,
                      addIndentLtx = FALSE, escPercLtx = FALSE, rndPer = 1,
                      rndNum1 = 2, rndNum2 = 1, msschar = "(missings)")
}

\arguments{
\item{x}{A data table. It will be turned into a data frame within the function using as.data.frame.}
\item{indx}{A list with which subgroups can be defined for which descriptive statistics should be generated. Per default, all observations are included. If indx has a length larger than one, then a merged descriptive table will be generated (see example below).}
\item{numDesc1}{The main descriptive statistic to be computed for numeric variables.}
\item{numDesc2}{The secondary descriptive statistic to be computed for numeric variables (presented in brackets).}
\item{showMiss}{Logical value indicating whether the numbers of missing observations should be presented.}
\item{showZeroMiss}{A locial value which determines whether zero missings should be indicated (per default missings are only shown if there are any). This option is only relevant when showMiss is TRUE.}
\item{boldVarnmsLtx}{A logical value which determines whether the variable names should be printed bold in a latex table.}
\item{addIndentLtx}{A logical value which determines whether there should be some indentation for the level names of factors in a latex table.}
\item{escPercLtx}{A logical value which determines whether for printing in latex the percent signs should be escaped using backward slashes.}
\item{rndPer}{The number of decimal places to round the relative proportions.}
\item{rndNum1}{The number of decimal places to round the main descriptive statistic of numeric variables.}
\item{rndNum2}{The number of decimal places to round the secondary descriptive statistic of numeric variables.}
\item{msschar}{The name of the table entry showing the numbers of missings.}
}


\details{
For non-numeric variables, the function automatically turns them into a
factor (using as.factor) and returns the number of cases and relative proportion
for each factor level. For numeric variables, the "quartile13" option returns
the first and third quartiles in brackets, while the "sd" option returns
the estimated standard deviation.


}

\value{
Returns a data frame containing the descriptive statistics.
}

\examples{
descTable(x = iris)
### For latex:
descTable(x = iris, boldVarnmsLtx = TRUE, addIndentLtx = TRUE,
          escPercLtx = TRUE)
### Multiple groups and showing missings:
a <- iris
a[100, 1] <- NA
descTable(x = a, indx = list("all obs" = 1:nrow(a),
          "only setosa" = which(a$Species=="setosa")))
}

