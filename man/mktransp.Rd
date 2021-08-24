\name{mktransp}
\alias{mktransp}
\title{Make a color transparent}
\description{
The function can be used to make a color value transparent with varying degrees of transparency. It can also be applied to vectors of colors.
}
\usage{
mktransp(color, alpha = 100)
}
\arguments{
\item{color}{Vector of colors to be turned transparent. Color values can be supplied in various forms.}

\item{alpha}{Level of transparency on a 0-255 scale. 255 refers to no transparency.}
}
\value{
The output is the transparent color(s) in hexadecimal notation.
}
\examples{
plot(1:4, pch=19, col=mktransp(c('blue', 'red')))
}

