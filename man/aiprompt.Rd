\name{aiprompt}
\alias{aiprompt}
\title{Send a prompt to a local AI model and get response}
\description{
The function sends a prompt to a local AI model and stores the answer
in a file.
}
\usage{
aiprompt(prmt, openFile = FALSE,
                     reqaddr="http://localhost:11434/api/generate",
                     modl="qwen3.5:9b",
                     temp_md = paste0(tempdir(), "/R_aipromptOutput.md"),
                     temp_pdf = paste0(tempdir(), "/R_aipromptOutput.pdf"),
                     timeoutsec = 120, ctxtwindow = 2^12)
}
\arguments{
\item{prmt}{The text prompt for the AI model.}
\item{openFile}{Logical, indicating whether the file containing the answer
should be automatically opened.}
\item{reqaddr}{The http address of the local AI model.}
\item{modl}{The name of the model in ollama to be used.}
\item{temp_md}{The name of the markdown file in which the generated
response is stored. If the file exists already, then the response is appended.}
\item{temp_pdf}{The name of the pdf file (rendered from the markdown file)
which contains the answer.}
\item{timeoutsec}{Seconds to wait before timeout happens.}
\item{ctxtwindow}{Amount of tokens defining the context size.}
}
\value{
Creates a file containing the answer which can be opened automatically.
}

