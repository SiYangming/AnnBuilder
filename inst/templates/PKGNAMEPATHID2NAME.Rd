\name{#PKGNAME#PATHID2NAME}
\alias{#PKGNAME#PATHID2NAME}
\title{An annotation data file that maps KEGG pathway identifiers to
  KEGG pathway names}
\description{
   #PKGNAME#PATHID2NAME maps KEGG pathway identifiers to pathway names used by
   KEGG for various pathways 
}
\details{
  This is an environment object containing key and value pairs. Keys are
  KEGG pathway identifiers and values are pathway names. Values are vectors of
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
	xx <- as.list(#PKGNAME#PATHID2NAME)
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

