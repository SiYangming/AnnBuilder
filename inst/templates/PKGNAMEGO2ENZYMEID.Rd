\name{#PKGNAME#GO2ENZYMEID}
\alias{#PKGNAME#GO2ENZYMEID}
\title{An annotation data file that maps Gene Ontology (GO) identifiers to
  Enzyme Commission numbers}
\description{
  #PKGNAME#GO2ENZYMEID maps GO ids to Enzyme Commission numbers
}
\details{
  This is an environment object containing key and value pairs. Keys are
  GO ids and values are Enzyme Commission numbers. Values are vectors of
  length 1. GO ids can not be mapped to any Enzyme Commission number are
  are assigned NAs.  

  Mappings are based on data provided by:

  #KEGGSOURCE#

  Package built: #DATE#
}
\references{
\url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#GO2ENZYMEID)
	if(length(xx) > 0){
		# Get the value of the first key
		xx[[1]]
		# Get values for a few keys
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

