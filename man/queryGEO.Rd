\name{queryGEO}
\alias{queryGEO}
\title{Function to extract a data file from the GEO web site}
\description{
  Data files that are available at GEO web site are identified by GEO
  accession numbers. Give a GEO object with the url for a common CGI and
  a GEO accession number, this function extracts data from the web site
  and returns a matrix containing the data portion of the file
}
\usage{
queryGEO(GEOObj, GEOAccNum)
}
\arguments{
  \item{GEOObj}{\code{GEOObj} a GEO object}
  \item{GEOAccNum}{\code{GEOAccNum} a character string for the GEO
    accession number of a desired file}
}
\details{
  The GEO object contains the url for a CGI script that processes user's
  request. \code{\link{queryGEO}} invokes the CGI by passing a GEO
  accession number and then processes the data file obtained.
}
\value{
  \code{\link{queryGEO}} returns a matrix containing data obtained.
}
\references{\url{www.ncbi.nlm.nih.gov/geo}}
\author{Jianhua Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{GEO-class}}}
\examples{
\dontrun{
geo <- GEO()
temp <- queryGEO(geo, "GPL49")
}
}
\keyword{manip}

