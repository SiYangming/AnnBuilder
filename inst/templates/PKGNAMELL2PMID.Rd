\name{#PKGNAME#LL2PMID}
\alias{#PKGNAME#LL2PMID}
\title{A data file that maps Entrez Gene identifiers to PubMed identifiers}
\description{
  #PKGNAME#LL2PMID maps Entrez Gene identifiers to PubMed identifiers representing
  publications that cited the genes corresponding to the Entrez Gene identifiers
}
\details{
  This is an environment object containing key and value pairs. Keys are
  Entrez Gene identifiers and values are the corresponding PubMed identifier. Values are
  vectors of length 1 or more depending on whether a Entrez Gene identifier can be
  mapped to one or more PubMed identifiers.

  Mappings were derived from data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#LL2PMID)
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

