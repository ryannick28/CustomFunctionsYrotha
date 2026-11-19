\name{renderoutp}
\alias{renderoutp}
\title{Render markdown text as a pdf file and open it.}
\description{
The function takes text ideally in markdown format and renders and opens it
as a pdf file.
}
\usage{
renderoutp(answ)
}
\arguments{
\item{answ}{A character element, ideally containing text in markdown format}
}


\examples{
# Don't run:
if(FALSE){
  res <- aiprompt(prmt = "Should you fertilize plants in winter")
  renderoutp(answ = res$response)
}
}
