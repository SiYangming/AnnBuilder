\name{UG-class}
\docType{class}
\alias{UG-class}
\alias{orgName}
\alias{orgName<-}
\alias{orgName,UG-method}
\alias{orgName<-,UG-method}
\alias{UG}
\title{Class "UG" a sub-class of pubRepo to handle data from UniGene}
\description{This class is a sub-class of pubRepo that is implemented
  specifically to parse data from UniGene (XX.data.gz, where XX is a
  abbreviation for a given organism)}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("UG", ...)}.
    A constructor (UG) is available and should be used to instatiate
    objects of this class
}
\section{Slots}{
  \describe{
    \item{\code{orgName}:}{Object of class \code{"character"} a
      character string for the name of the organism of concern}
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"pubRepo"} a character string for the url of the source data}
    \item{\code{parser}:}{Object of class \code{"character", from class
	"pubRepo"}  a character string for the name of the file
      containing a segment of perl code with instructions on how the
      source data will be processed and output be generated} 
    \item{\code{baseFile}:}{Object of class \code{"character", from
	class "pubRepo"} a character string for the name of the file
      that contains data that will be used as the base to process the
      source data.  Data from the source that are related to elements in
      the base file will be extracted. baseFile is assumed to be a two
      folumn file with the first column being some type of arbitrary ids
      (e.g. Affymetrix probe ids) and the second cloumn being the
      corresponding ids of a given public repository (e.g. GenBank
      accession numbers or UniGene ids)}
  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
  \describe{
    \item{orgName<-}{\code{signature(object = "UG")}: Sets the value
      for the organism slot}
    \item{orgName}{\code{signature(object = "UG")}: Gets the value for
      the organism slot}
  }
}
\references{\url{www.ncbi.nlm.nih.gov/UniGene}}
\author{Jianhua Zhang}

\seealso{\code{\link{pubRepo-class}}}

\examples{
\dontrun{
# Parse a truncated version of Hs.data.gz from Bioconductor
path <- file.path(.path.package("pubRepo"), "data")
temp <- matrix(c("32469_f_at", "D90278", "32469_at", "L00693", "33825_at",
"X68733", "35730_at", "X03350", "38912_at", "D90042", "38936_at",
"M16652"), ncol = 2, byrow = TRUE)
write.table(temp, "tempfile", sep = "\t", quote = FALSE,
row.names = FALSE, col.names = FALSE)  
ug <- UG(srcUrl =
"http://www.bioconductor.org/datafiles/wwwsources/Ths.data.gz",
parser = file.path(path, "basedUGParser"), baseFile = "tempfile",
organism = "human")
data <- parseData(ug)
unlink("tempfile")
}
}
\keyword{classes}
