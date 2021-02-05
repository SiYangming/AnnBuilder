\name{#PKGNAME#PATHID2EXTID}
\alias{#PKGNAME#PATHID2EXTID}
\title{An annotation data file that maps KEGG pathway ids to LocusLink
  or Open Reading Frame ids}
\description{
  #PKGNAME#PATHID2EXTID maps KEGG pathway ids to LocusLink (for human,
  mouse, and rat) or Open Reading Frame (for yeast) ids   
}
\details{
  This is an environment object containing key and value pairs. Keys are
  KEGG pathway ids and values are LocusLink ids for human, mouse, and
  rat or Open Reading Frame (ORF) ids for yeast. Values are
  vectors of length 1 or greater depending on whether a pathway id can
  be maapped to one or more LocusLink or ORF ids. NAs are
  assigned to KEGG pathway ids that can not be mapped to any LocusLink
  or ORF ids. 

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
	xx <- as.list(#PKGNAME#PATHID2EXTID)
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

