\name{GOCCOFFSPRING}
\alias{GOCCOFFSPRING}
\title{Annotation of GO Identifiers to their Cellular Component Offspring} 
\description{
This is an R environment (hash table) mapping GO ids to all the GO
ids of their offspring based on the  cellular component (CC) GO terms. Offspring are
defined as the direct or indirect children/grand children nodes of a
given node for a GO id based on the directed acyclic graph defined by
Gene Ontology Consortium. GO terms for offspring nodes are more specific
than that of the parent.
}
\details{
  GO ids are keys and the corresponding offspring GO ids are
  values. Values are vectors of GO ids for all the offspring of the key
  GO id. GO ids that do not have any offspring node are assigned NA as
  the value. 

  Cellular component is defined as the subcellular structures,
  locations, and macromolecular complexes; examples include nucleus,
  telomere, and origin recognition complex as defined b y Gene Ontology
  Consortium. 
  
  Mappings were based on data provided:

  #GOSOURCE#
 
  Package built: #DATE#
}
\references{
\url{http://www.geneontology.org/} and
\url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} 
}
\examples{
  require("GO", character.only = TRUE) || stop("GO unavailable")
  # Convert the environment object to a list
  xx <- as.list(GOCCOFFSPRING)
  # Remove GO ids that do not have any offspring
  xx <- xx[!is.na(xx)]
  if(length(xx) > 0){
    # Get the offspring GO ids for the first two elents of xx
    goids <- xx[1:2]
  }
}
\keyword{datasets}

