\name{#PKGNAME#DESCRIPTION}
\alias{#PKGNAME#DESCRIPTION}
\title{An annotation data file that maps Open Reading Frame (ORF)
  identifiers to textural descriptions of the corresponding genes}
\description{
  #PKGNAME#DESCRIPTION maps yeast ORF ids to descriptive information
  about genes corresponding to the ORF ids 
}
\details{
  This is an environment object containing key and value pairs. Keys are
  ORF ids and values are the corresponding descriptions of genes. Values
  are vectors of length 1. Probe ids that can not be mapped to
  descriptive information are assigned a value NA.
  
  Annotation based on data provided by:

  #YGSOURCE#

  Package built: #DATE#
}
\references{
\url{http://www.yeastgenome.org/DownloadContents.shtml}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#DESCRIPTION)
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

