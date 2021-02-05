\name{IPI-class}   
\docType{class}
\alias{IPI-class}
\alias{parseData,IPI-method}
\alias{IPI}
\title{Class "IPI" a sub-class of pubRepo to handle data from
  International Protein Index (IPI)}
\description{This class is a sub-class of pubRepo that is implemented
  specifically to parse data from IPI (ipi.*.dat.gz)}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IPI", ...)}.
    A constructor (IPI) is available and should be used to instatiate
    objects of \code{\link{IPI}} 
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
	class "pubRepo"} a character string for the name of the gzipped file
      that contains data from IPI ftp site.  For example,
      ipi.HUMAN.dat.gz is the file for human, and ipi.MOUSE.dat.gz is
      for mouse, etc.}
  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
    \describe{
    \item{parseData}{\code{signature(object = "IPI")}: A method to parse
      a source file using a specified parser}
  }

}
\references{\url{http://www.ebi.ac.uk/IPI/IPIhelp.html}}
\author{Ting-Yuan Liu}

\seealso{\code{\link{pubRepo-class}}}

\examples{
\dontrun{
## create IPI class
ipi <- IPI(srcUrl="ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/",
           organism = "human")

## Parse ipi.HUMAN.dat.gz from IPI ftp site
tmpFile <- loadFromUrl(paste(srcUrls(ipi),baseFile(ipi)), sep="")
system("grep "//" ")
con <- file(tmpFile, "r")

tmpRead <- readLines(con, n=200)
endSymbol <- grep("//", tmpRead)
tmpRead <- tmpRead[1:endSymbol[length(endSymbol)]]
file <- tempfile()
writeLines(tmpRead, file)

system(paste("mv ", file, " ", tempdir(), "/ipi.tiny.dat", sep=""))
system(paste("gzip ", tempdir(), "/ipi.tiny.dat", sep=""))
ipiParser(ipiData=paste(tempdir(), "/ipi.tiny.dat.gz", sep=""), fromWeb=FALSE)

}
}
\keyword{classes}
