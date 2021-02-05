\name{#PKGNAME#PATHNAME2ID}
\alias{#PKGNAME#PATHNAME2ID}
\title{Map KEGG pathway names to pathway ids}
\description{
  #PKGNAME#PATHNAME2ID is an R environment that maps KEGG pathway names
  to the ids used by KEGG for various pathways 
}
\details{
  Each KEGG pathway name is mapped to a vector containing a single id
  for the KEGG pathway.
  
  Mappings were based on data provided by:

  #KEGGSOURCE#
  
  Package built: #DATE#
}
\references{
\url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#PATHNAME2ID)
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

