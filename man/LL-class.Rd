\name{LL-class}   
\docType{class}
\alias{LL-class}
\alias{LL}
\title{Class "LL" a sub-class of pubRepo to handle data from LocusLink}
\description{This class is a sub-class of pubRepo that is implemented
  specifically to parse data from LocusLink (ll\_teml.gz)}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LL", ...)}.
    A constructor (LL) is available and should be used to instatiate
    objects of \code{\link{LL}} 
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
No methods defined with class "LL" in the signature.
}
\references{\url{www.ncbi.nlm.nih.gov/LocusLink}}
\author{Jianhua Zhang}
\note{This class is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics funtionalities through R}

\seealso{\code{\link{pubRepo-class}}}

\examples{
\dontrun{
# Parse a truncated version of LL\_tmpl.gz from Bioconductor
path <- file.path(.path.package("AnnBuilder"), "scripts")
temp <- matrix(c("32469_f_at", "D90278", "32469_at", "L00693", "33825_at",
"X68733", "35730_at", "X03350", "38912_at", "D90042", "38936_at",
"M16652"), ncol = 2, byrow = TRUE)
write.table(temp, "tempfile", sep = "\t", quote = FALSE,
row.names = FALSE, col.names = FALSE)  
ll <- LL(srcUrl =
"http://www.bioconductor.org/datafiles/wwwsources/Tll_tmpl.gz",
parser = file.path(path, "gbLLParser"), baseFile = "tempfile")
data <- parseData(ll)
unlink("tempfile")
}
}
\keyword{classes}
