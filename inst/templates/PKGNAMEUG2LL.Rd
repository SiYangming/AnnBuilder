\name{#PKGNAME#UG2LL}
\alias{#PKGNAME#UG2LL}
\title{An annotation data file that maps UniGene cluster identifiers to
  Entrez Gene identifiers} 
\description{
  #PKGNAME#UG2LL maps UniGene identifiers to Entrez Gene identifiers
}
\details{
  This is an environment object containing key and value pairs. Keys are
  UniGene identifiers and values are Entrez Gene identifiers. Values are vectors of length
  1.

  Mappings were derived from data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#UG2LL)
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

