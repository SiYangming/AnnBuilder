\name{#PKGNAME#GO2ALLPROBES}
\alias{#PKGNAME#GO2ALLPROBES}
\title{Map Between Gene Ontology (GO) Identifiers and all Manufacturer
  Identifiers in the subtree}
\description{
  #PKGNAME#GO2ALLPROBES is an R environment that provides mappings
  between a given GO identifier and all manufactuerer identifiers annotated at that
  GO term or one of its children in the GO ontology.
}
\details{
  GO consists of three ontologies---molecular function (MF),
  biological process (BP), and cellular component (CC).  All ontologies
  are structured as directed acyclic graphs (DAGs).  Each node in
  each DAG (tree) is a GO term (id) associated with a named vector of
  manufacturer identifiers.  The name associated with each manufacturer id
  corresponds to the evidence code for that GO identifier.  This environment
  #PKGNAME#GO2ALLPROBES maps between a given GO identifier and all manufactuerer
  identifiers annotated at that GO term or one of its children in the GO ontology.  

  The evidence code indicates what kind of evidence supports the
  association between the GO and Entrez Gene identifiers. Evidence codes
  currently in use include:

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

  A GO identifier may be mapped to the same manufacturer identifier more than once but the
  evidence code can be different. Mappings between Gene Ontology identifiers and
  Gene Ontology terms and other information are available in a separate
  data package named GO.  
  
  Mappings were based on data provided by:

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
        # Gets the probe identifiers for the top 2nd and 3nd GO identifiers
        goids <- xx[2:3]
        # Gets all the probe identifiers for the first element of goids
        goids[[1]]
        # Evidence code for the mappings
        names(goids[[1]])
    }
}
\keyword{datasets}

