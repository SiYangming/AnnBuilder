\name{YEAST-class}   
\docType{class}
\alias{YEAST-class}
\alias{parseData,YEAST-method}
\alias{YEAST}
\title{Class "YEAST" a sub-class of pubRepo to handle data from ftp.yeastgenome.org}
\description{This class is a sub-class of pubRepo that is implemented
  specifically to parse the data
  \url{ftp://ftp.yeastgenome.org/pub/yeast/sequence_similarity/domains/domains.tab}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("YEAST", ...)}.
    A constructor (YEAST) is available and should be used to instatiate
    objects of \code{\link{YEAST}} 
}
\section{Slots}{
  \describe{
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"pubRepo"} a character string for the surce url where data
      will be downloaded/processed}
    \item{\code{parser}:}{Object of class \code{"character", from class
	"pubRepo"}  a character string for the name of the file
      containing a segment of perl code with instructions on how the
      source data will be processed and output be generated}
    \item{\code{baseFile}:}{Object of class \code{"character", from
	class "pubRepo"} a character string for the name of the file
      used to be parsed.  The default file name is "domains.tab".  }
  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
    \describe{
    \item{parseData}{\code{signature(object = "YEAST")}: A method to parse
      a source file using a specified parser}
  }

}
\references{\url{ftp://ftp.yeastgenome.org/pub/yeast/sequence_similarity/domains/domains.tab}}
\author{Ting-Yuan Liu}
\seealso{\code{\link{pubRepo-class}}}

\examples{
\dontrun{
## create YEAST class
yeast <- YEAST(srcUrl="ftp://ftp.yeastgenome.org/pub/yeast/sequence_similarity/domains/",
           baseFile="domains.tab")
yeastDomain <- parseData(yeast)
head(yeastDomain)
}
}
\keyword{classes}
