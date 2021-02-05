\name{#PKGNAME#PATHNAME2ID} 
\alias{#PKGNAME#PATHNAME2ID}
\title{An annotation data file that maps KEGG pathway names to
  identifiers for the corresponding pathway names used by KEGG}
\description{
   #PKGNAME#PATHNAME2ID maps KEGG pathway names to pathway identifiers used by
   KEGG for various pathways 
}
\details{
  This is an environment object containing key and value pairs. Keys are
  KEGG pathway names and values are pathway identifiers. Values are vectors of
  length 1.  

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
		# get the value for the first key
		xx[[1]]
		# Get the values for a few keys
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

