\name{#PKGNAME#EXTID2PATHID}
\alias{#PKGNAME#EXTID2PATHID}
\title{An annotation data file that maps Entrez Gene or Open Reading Frame
  identifiers KEGG pathway identifiers}
\description{
  #PKGNAME#EXTID2PATHID maps Entrez Gene (for human, mouse, and rat)
  and Open Reading Frame (for yeast) identifiers to KEGG pathway identifiers.   
}
\details{
  This is an environment object containing key and value pairs. Keys are
  Entrez Gene identifiers for human, mouse, and rat and Open Reading Frame (ORF) identifiers for
  yeast and values are the corresponding KEGG pathway identifiers. Values are
  vectors of length 1 or greater depending on whether a given external
  identifier can be mapped to only one or more KEGG pathway identifiers.  NAs are
  assigned to Entrez Gene or ORF identifiers that can not be mapped to any pathway
  identifier.

  KEGG pathway identifiers are the identifiers used by KEGG for various
  pathways.

  Mappings between KEGG pathway identifiers and pathway names can be obtained
  through another environment named PATHID2NAME. 

  Mappings were based on data provided by:

  #KEGGSOURCE#

  Package built: #DATE#
}
\references{
\url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#EXTID2PATHID)
	if(length(xx) > 0){
		# Get the value of the first key
		xx[[1]]
		# Get the values for multiget for a few keys
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

