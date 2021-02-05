\name{PFAM-class}   
\docType{class}
\alias{PFAM-class}
\alias{parseData,PFAM-method}
\alias{PFAM}
\title{Class "PFAM" a sub-class of pubRepo to handle data from http://www.sanger.ac.uk/Software/Pfam/}
\description{This class is a sub-class of pubRepo that is implemented
  specifically to parse the data
  \url{ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PFAM", ...)}.
    A constructor (PFAM) is available and should be used to instatiate
    objects of \code{\link{PFAM}} 
}
\section{Slots}{
  \describe{
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"pubRepo"} a character string for the surce url where data
      will be downloaded/processed}
    \item{\code{fromWeb}:}{Object of class \code{"logical"} for
      inficating whether the source file should be accessed through the
      web or locally}

  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
    \describe{
    \item{parseData}{\code{signature(object = "PFAM")}: A method to parse
      a source file using a specified parser}
  }

}
\references{\url{ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz}}
\author{Ting-Yuan Liu}
\seealso{\code{\link{pubRepo-class}}}

\examples{
\dontrun{
## create PFAM class
pfamObj <- PFAM(srcUrl="ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz",
                fromWeb=TRUE)
tableList <- parseData(pfamObj)
}
}
\keyword{classes}
