\name{#PKGNAME#GO2LL}
\alias{#PKGNAME#GO2LL}
\title{An annotation data that maps Gene Ontology identifiers to
  LocusLink identifiers} 
\description{
  #PKGNAME#GO2LL maps Gene Ontology (GO) ids to LocusLink ids for genes
  whose products have the functions corresponding to the GO ids
}
\details{
  This is an environment object containing key and value pairs. Keys are
  GO ids and values are LocusLink ids. Values are vectors of length 1 or
  more depending on whether a GO id can be mapped to one or more
  LocusLink ids. Names of values are evidence code that indicate the
  type of evidences that support associations between genes and GO
  functions. 

  The evidence codes in use include:
  
  IMP: inferred from mutant phenotype  
  IGI: inferred from genetic interaction
  IPI: inferred from physical interaction  
  ISS: inferred from sequence similarity 
  IDA: inferred from direct assay  
  IEP: inferred from expression pattern  
  IEA: inferred from electronic annotation  
  TAS: traceable author statement  
  NAS: non-traceable author statement  
  ND: no biological data available  
  IC: inferred by curator

  Mappings were derived from data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{ 
\url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#GO2LL)
	if(length(xx) > 0){
		# Get the value of the first key
		xx[[1]]
		# Get the values  for a few keys
		if(length(xx) >= 3){
			xx[1:3]
		}
	}
}
\keyword{datasets}

