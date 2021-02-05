\name{#PKGNAME#PMID2LL}
\alias{#PKGNAME#PMID2LL}
\title{An annotation data file that maps PubMed identifiers to LocusLink
  identifiers} 
\description{
  #PKGNAME#PMID2LL maps PubMed ids to LocusLink ids for genes that have
  been cited by publications represented by the PubMed ids
}
\details{
  This is an environment object containing key and value pairs. Keys are
  PubMed id and values are LocusLink ids representing genes cited by
  publications corresponding to the PubMed ids. Values are vectors of
  length 1 or more depending on whether a PubMed id can be mapped to one
  or more LocusLink ids.

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

