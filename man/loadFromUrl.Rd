\name{loadFromUrl}
\alias{loadFromUrl}
\alias{validateUrl}
\alias{unzipFile}
\title{Functions to load files from a web site}
\description{
  Given an url, these functions download a file from a given web site
  and unzip the file if it is compressed.
}
\usage{
loadFromUrl(srcUrl, destDir = "", verbose=FALSE)
validateUrl(srcUrl)
unzipFile(fileName, where = file.path(.path.package("AnnBuilder"),
"data"), isgz = FALSE)
}

\arguments{
  \item{srcUrl}{\code{srcUrl} a character string for the url of the file
    to be downloaded}
  \item{destDir}{\code{destDir} a character string for a loacal
    directory where the file to be downloaded will be saved}
  \item{where}{\code{where} same as destDir}
  \item{isgz}{\code{isga} a boolean indicating whether the downloaded
    file is a gz file}
  \item{fileName}{\code{fileName} a character string for the name of a
    file} 
  \item{verbose}{A booline indicating whether to print extra information.}
}
\details{
  These functions are used by various objects in package pubRepo to
  download data files from a web site. If the file is compressed,
  decompressing will be applied and the path for the decompressed file
  will be returned.

  \code{\link{validateUrl}} will terminate the process if an invalid url
  is passed.

  \code{\link{unzipFile}} decompress the file passed as fileName.
}
\value{
  \code{\link{loadFromUrl}} returns a character string for the name of
  the file saved locally.
}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{pubRepo-class}}}
\examples{
\dontrun{
# Get a dummy data file from Bioconductor web site
data <-
loadFromUrl("http://www.bioconductor.org/datafiles/wwwsources/Tll_tmpl.gz",
destDir = "")
unlink(data)
}
}
\keyword{manip}

