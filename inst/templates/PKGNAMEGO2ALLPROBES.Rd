\name{#PKGNAME#GO2ALLPROBES}
\alias{#PKGNAME#GO2ALLPROBES}
\title{Map Between Gene Ontology (GO) Identifiers and all Manufacturer Identifiers}
\description{
  #PKGNAME#GO2ALLPROBES is an R environment that provides mappings
  between a given GO id and all manufacturer ids.
}
\details{
  Each GO term maps to a named vector of manufacturer ids. The name
  associated with each manufacturer id corresponds to the evidence code
  for that GO id. The evidence code indicates what kind of evidence
  supports the association between the GO and Entrez Gene ids. Evidence
  codes currently in use include:

  IMP - inferred from mutant phenotype  

  IGI - inferred from genetic interaction

  IPI - inferred from physical interaction  

  ISS - inferred from sequence similarity  

  IDA - inferred from direct assay  

  IEP - inferred from expression pattern  

  IEA - inferred from electronic annotation  

  TAS - traceable author statement  

  NAS - non-traceable author statement  

  ND - no biological data available  

  IC - inferred by curator

  A GO id may be mapped to the same manufacturer id more than once but the
  evidence code can be different. Mappings between Gene Ontology ids and
  Gene Ontology terms and other information are available in a separate
  data package named GO.  
  
  Mappings were based on data provided:

  #LLSOURCE#
  
  Package built #DATE#
}
\references{
\url{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/}
}
\examples{
    # Convert to a list
    xx <- as.list(#PKGNAME#GO2ALLPROBES)
    if(length(xx) > 0){
        # Gets the probe ids for the top 2nd and 3nd GO ids
        goids <- xx[2:3]
        # Gets all the probe ids for the first element of goids
        goids[[1]]
        # Evidence code for the mappings
        names(goids[[1]])
    }
}
\keyword{datasets}

