\name{#PKGNAME#ACCNUM2LL}
\alias{#PKGNAME#ACCNUM2LL}
\title{An annotation data file that maps GenBank Accession numbers to
  Entrez Gene identifiers}
\description{
  #PKGNAME#ACCNUM2LL maps GenBank Accession numbers to Entrez Gene identifiers
}
\details{
  This is an environment object containing key and value pairs. Keys are
  GenBank Accession numbers and values are the corresponding Entrez Gene
  identifiers. Values are vectors of length 1 or more depending on whether a
  GenBank Accession number can be mapped to one or more Entrez Gene identifiers.

  Mappings were derived from data provided by:

  #LLSOURCE#

  #UGSOURCE# 
  
  Package built: #DATE#
}
\references{
\url{http://www.ncbi.nlm.nih.gov}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#ACCNUM2LL)
	if(length(xx) > 0){
		# Get the value of the first key
		xx[1]
		# Get the values for a few keys
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

