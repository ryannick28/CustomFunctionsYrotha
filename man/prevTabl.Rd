\name{prevTabl}
\alias{prevTabl}
\title{Create prevalence table}
\description{
This function creates a prevalence table for a tested condition in a data table. One can define a function which returns a logic output for each element in the data. The function then counts how many cases have returned TRUE for each column and presents this information. In addition, one can define a categorical variable and return the results for each level separately. The function can deal with missing values in the data.}
\usage{
prevTabl(X, FUN, catVar = NULL, atLeastOnce = FALSE)
}
\arguments{
\item{X}{A data table, classically a data frame.}
\item{FUN}{A function to test a specific condition along each column of the data. Must return a logical value for each element (or NA).}
\item{catVar}{A factor with the same length as there are rows in the data. When supplied, the prevalences of successes will be presented for each level of catVar separately.}
\item{atLeastOnce}{A logical value. If set to TRUE, the final table will include an additional column, which will present whether a specific case scored at least one TRUE over all columns of the data.}
}
\value{
The (silent) return object comprises the number of successes, the number of observations and the successrate for each column. If atLeastOnce was set to TRUE, then the same information is additionally returned for this extra column. The final prevalence table is printed when the function is called and also stored in the return object.
}
\examples{
f <- function(x){x < 3.1}
res <- prevTabl(X = iris[,-5], FUN = f, catVar = iris$Species, atLeastOnce = TRUE)
res
}

