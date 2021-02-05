\name{GEO-class}
\docType{class}
\alias{GEO-class}
\alias{GEO}
\alias{readData,GEO-method}
\title{Class "GEO" represents a GEO object that reads/downloads data
  from the GEO web site}
\description{The GEO web site contains data files represented by GEO
  accession numbers. Class GEO reads/downloads data files from the site
  if correct url and GEO accession numbers are provided}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("UG", ...)}.
    A constructor (GEO) is available and should be used to instantiate
    objects of this class  
}
\section{Slots}{
  \describe{
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"pubRepo"} - a character string for the url of a CGI script that
    handles data requests, which is:
    \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?} at the time of
    writing} 
  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
  \describe{
    \item{readData}{\code{signature(object = "GEO")}: reads data from
      GEO and then parses the data to a matrix}
  }
}
\references{Programming with data}
\author{Jianhua Zhang}

\seealso{\code{\link{queryGEO}},\code{\link{pubRepo-class}}}

\examples{
\dontrun{
geo <- GEO("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?")
# The GEOAccNum may be invalid due to changes at GEO site
data <- readData(geo, GEOAccNum = "GPL16" )
}}
\keyword{classes}


