\name{mmdgp}
\alias{mmdgp}
\title{Generate longitudinal hierarchical data}
\description{
A function to generate longitudinal hierarchical data. The data represents measurements from people, nested in classes. Each person is measured once per day for a number of days. There is an age, a gender and a day effect. The random effects are only random intercepts for people and class. Can for example be used to test some analysis or some plotting procedure.
}
\usage{
mmdgp(n = 200, nC = 5, sd_S = 3, sd_C = 5, sd_e = 2,
b0 = 50, tb = c(0,4,10,4,0), genb = 10, ageb = 1)
}
\arguments{
\item{n}{Number of people}
\item{nC}{Number of classes}
\item{sd_S}{Standard deviation of people (random intercept)}
\item{sd_C}{Standard deviation of classes (random intercept)}
\item{sd_e}{Standard deviation of errors}
\item{b0}{Intercept}
\item{tb}{A vector containing the effects of each day. Length can be chosen depending on how many days should be included in the data set.}
\item{genb}{Gender effect}
\item{ageb}{Age effect}
}
\value{
A data frame with the predictor variables and the target variable.
}
\examples{
dat = mmdgp()
# Lmer fit:
# lmerTest::lmer(y ~ gen + age + day + (1|id) + (1|class), dat)
}

