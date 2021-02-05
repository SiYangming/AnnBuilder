\name{#PKGNAME#GO}
\alias{#PKGNAME#GO}
\title{Map between Manufacturer IDs and Gene Ontology (GO)}
\description{
  #PKGNAME#GO is an R environment that provides mappings between
  manufacturer ids and GO ids.
}
\details{
  Each manufacturer id is mapped to a named vector of GO ids. The name
  associated with each GO id corresponds to the evidence code for that
  GO id. Evidence codes currently in use include:
  
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
  
  Mappings between probe ids and GO information were obtained through
  their mappings to Entrez Gene ids. NAs are assigned to probe identifiers
  that can not be mapped to any Gene Ontology information. Mappings
  between Gene Ontology ids an Gene Ontology terms and other information
  are available in a separate data package named GO.
  
  Mappings were based on data provided by:
  
  #LLSOURCE#

  #GOSOURCE#
  
  Package built #DATE#
}

\examples{
    # Convert to a list
    xx <- as.list(#PKGNAME#GO)
    # Remove all the NAs
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
            # Try the firest one
            got <- xx[[1]]           
            got[[1]][["GOID"]]
            got[[1]][["Ontology"]]
            got[[1]][["Evidence"]]

    }
}
\keyword{datasets}

