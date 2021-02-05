\name{#PKGNAME#ENZYMEID2GO}
\alias{#PKGNAME#ENZYMEID2GO}
\title{An annotation data file that maps Enzyme Commission numbers to
  Gene Ontology identifiers}
\description{
  #PKGNAME#ENZYMEID2GO maps Enzyme Commission numbers to Gene Ontoloty
  (GO) identifiers 
}
\details{
  This is an environment object containing key and value pairs. Keys are
  Enzyme Commission numbers and values are GO identifiers. Values are vectors of
  length 1. Enzyme Commission numbers that can not be mapped to a GO
  term are assigned a value NA.

  Mappings were based on data provided by:

  #KEGGSOURCE#
  
  Package built: #DATE#
}
\references{
\url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#ENZYMEID2GO)
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

