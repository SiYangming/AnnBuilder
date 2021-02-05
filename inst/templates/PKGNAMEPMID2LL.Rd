\name{#PKGNAME#PMID2LL}
\alias{#PKGNAME#PMID2LL}
\title{An annotation data file that maps PubMed identifiers to Entrez Gene
  identifiers} 
\description{
  #PKGNAME#PMID2LL maps PubMed identifiers to Entrez Gene identifiers for genes that have
  been cited by publications represented by the PubMed identifiers
}
\details{
  This is an environment object containing key and value pairs. Keys are
  PubMed identifier and values are Entrez Gene identifiers representing genes cited by
  publications corresponding to the PubMed identifiers. Values are vectors of
  length 1 or more depending on whether a PubMed identifier can be mapped to one
  or more Entrez Gene identifiers.

  Mappings were derived from data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#PMID2LL)
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

