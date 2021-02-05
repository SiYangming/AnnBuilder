\name{#PKGNAME#ENTREZID2CHR}
\alias{#PKGNAME#ENTREZID2CHR}
\alias{#PKGNAME#LOCUSID2CHR}
\title{An annotation data file that maps Entrez Gene identifiers to
  chromosome number}
\description{
  #PKGNAME#ENTREZID2CHR maps Entrez Gene identifiers to the chromosome numbers the
  genes represented by the Locuslink identifiers reside
}
\details{
  This is an environment object containing key and value pairs. Keys are
  Entrez Gene identifiers and values are the corresponding chromosome numbers the
  genes reside. Values are vectors of length 1 or more depending on
  whether a give Entrez Gene identifier can be mapped to one or more chromosomes.

  Mappings were derived from data provided by:

  #GPSOURCE#

  Package built: #DATE#
}
\references{
\url{http://www.genome.ucsc.edu/goldenPath/hg16/database/}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#ENTREZID2CHR)
	if(length(xx) > 0){
		# Get the value of the first Entrez Gene id
		xx[1]
		# Get the values for a few Entrez Gene identifiers
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

