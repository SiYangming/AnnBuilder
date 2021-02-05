\name{#PKGNAME#LL2UG}
\alias{#PKGNAME#LL2UG}
\title{An annotation data file that maps LocusLink identifiers to
  UniGene cluster identifiers}
\description{
  #PKGNAME#LL2UG maps LocusLink ids to UniGene ids
}
\details{
  This is an environment object containing key and value pairs. Keys are
  LocusLink ids and values are the corresponding UniGene ids. Values are
  vectors of length 1 or more depending on whether a LocusLink id can be
  mapped to one or more UniGene ids.

  Mappings were derived from data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#LL2UG)
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

