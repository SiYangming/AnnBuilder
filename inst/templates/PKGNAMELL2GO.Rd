\name{#PKGNAME#LL2GO}
\alias{#PKGNAME#LL2GO}
\title{An annotation data file that maps Entrez Gene identifiers to Gene
  Ontology identifiers} 
\description{
  #PKGNAME#LL2GO maps Entrez Gene identifiers to Gene Ontology (GO) identifiers
  representing functions of the products of genes correspinding to the
  Entrez Gene identifiers 
}
\details{
  This is an environment object containing key and value pairs. Keys are
  Entrez Gene identifiers and values are GO identifiers. Values are vectors of length 1 or
  more depending on whether a Entrez Gene identifier can be mapped to one or more
  Gene Ontology identifiers. Names of values are evidence code.

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
  ND: no biological data avai

  Mappings were derived from data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME#LL2GO)
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

