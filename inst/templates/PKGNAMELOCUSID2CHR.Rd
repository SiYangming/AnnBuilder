\name{#PKGNAME#LOCUSID2CHR}
\alias{#PKGNAME#LOCUSID2CHR}
\title{An annotation data file that maps LocusLink identifiers to
  chromosome number}
\description{
  #PKGNAME#LOCUSID2CHR maps LocusLink ids to the chromosome numbers the
  genes represented by the Locuslink ids reside
}
\details{
  This is an environment object containing key and value pairs. Keys are
  LocusLink ids and values are the corresponding chromosome numbers the
  genes reside. Values are vectors of length 1 or more depending on
  whether a give LocusLink id can be mapped to one or more chromosomes.

  Mappings were derived from data provided by:

  #GPSOURCE#

  Package built: #DATE#
}
\references{
\url{http://www.genome.ucsc.edu/goldenPath/hg16/database/}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#LOCUSID2CHR)
	if(length(xx) > 0){
		# Get the value of the first LocusLink id
		xx[1]
		# Get the values for a few LocusLink ids
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

