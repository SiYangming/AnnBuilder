\name{#PKGNAME#ORF}
\alias{#PKGNAME#ORF}
\title{Map Manufacturer Identifiers to Open Reading Frame (ORF) Identifiers} 
\description{
  #PKGNAME#ORF is an R environment that provides mappings between
  manufacturer and ORF identifiers.
}
\details{
  Each manufacturer identifier is mapped to a vector of ORF identifiers. The length of
  the vector may be one or longer, depending on how many ORF identifiers the
  manufacturer identifier can be mapped to. An \code{NA} is reported for any
  manufacturer identifier that cannot be mapped to an ORF identifier at this time.
  
  Annotation based on data provided by:

  #YGSOURCE#
  
  Package built: #DATE#
}
\references{
\url{http://www.yeastgenome.org/DownloadContents.shtml}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#ORF)
	if(length(xx) > 0){
		# Get the value of the first key
		xx[[1]]
		# Get the values for a few keys
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

