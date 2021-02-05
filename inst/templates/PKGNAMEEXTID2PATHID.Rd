\name{#PKGNAME#EXTID2PATHID}
\alias{#PKGNAME#EXTID2PATHID}
\title{An annotation data file that maps LocusLink or Open Reading Frame
  ids KEGG pathway ids}
\description{
  #PKGNAME#EXTID2PATHID maps LocusLink (for human, mouse, and rat)
  and Open Reading Frame (for yeast) ids to KEGG pathway ids.   
}
\details{
  This is an environment object containing key and value pairs. Keys are
  LocusLink ids for human, mouse, and rat and Open Reading Frame (ORF) ids for
  yeast and values are the corresponding KEGG pathway ids. Values are
  vectors of length 1 or greater depending on whether a given external
  id can be mapped to only one or more KEGG pathway ids.  NAs are
  assigned to LocusLink or ORF ids that can not be mapped to any pathway
  id.

  KEGG pathway ids are the identifiers used by KEGG for various
  pathways.

  Mappings between KEGG pathway ids and pathway names can be obtained
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

