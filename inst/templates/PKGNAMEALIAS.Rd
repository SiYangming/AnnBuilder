\name{#PKGNAME#ALIAS}
\alias{#PKGNAME#ALIAS}
\title{Map Open Reading Frame (ORF) Identifiers to Alias Gene Names}
\description{
  A set of gene names may have been used to report yeast genes
  represented by ORF identifiers. One of these names is chosen to be the primary
  gene name, and the others are considered aliases. This R environment
  provides mappings between the primary name and aliases.
}
\details{
  Each primary name maps to a vector of alias names. If there are no
  aliases, the vector will contain \code{NA}.
  
  Annotation based on data provided by:

  #YGSOURCE#
  
  Package built: #DATE#
}
\references{
\url{http://www.yeastgenome.org/DownloadContents.shtml}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#ALIAS)
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

