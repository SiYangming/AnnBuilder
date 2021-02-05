\name{#PKGNAME#LL2ACCNUM}
\alias{#PKGNAME#LL2ACCNUM}
\title{An annotation data file that maps Entrez Gene identifiers to
  GenBank Accession numbers}
\description{
  #PKGNAME#LL2ACCNUM maps Entrez Gene identifiers to GenBank Accession numbers
}
\details{
  This is an environment object containing key and value pairs. Keys are
  Entrez Gene identifiers and values are the corresponding GenBank Accession
  Numbers. Values are vectors of length 1 or more depending on whether a
  Entrez Gene identifier can be mapped to one or more GenBank Accession Numbers.

  Mappings were derived from data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{http://www.ncbi.nlm.nih.gov}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#LL2ACCNUM)
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

