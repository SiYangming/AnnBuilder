\name{GOCCCHILDREN}
\alias{GOCCCHILDREN}
\title{Annotation of GO Identifiers to their Cellular Component Children} 
\description{
  This data set describes associations between GO molecular function (CC)
  terms and their direct children CC terms, based on the directed acyclic
  graph (DAG) defined by the Gene Ontology Consortium. The format is an R
  environment mapping the GO CC terms to all direct children terms, where a
  direct child term is a more specific GO term that is immediately preceded
  by the given GO term in the DAG.
}
\details{
  Each GO CC term is mapped to a vector of children GO CC terms.

  Cellular component is defined as the subcellular structures,
  locations, and macromolecular complexes; examples include nucleus,
  telomere, and origin recognition complex as defined by Gene Ontology
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
  require("annotate", character.only = TRUE) || stop("annotate unavailable")
  # Convert the environment object to a list
  xx <- as.list(GOCCCHILDREN)
  # Remove GO IDs that do not have any children
  xx <- xx[!is.na(xx)]
  if(length(xx) > 0){
    goids <- xx[[1]]
    # Find out the GO terms for the first parent goid
    GOID(GOTERM[[goids[1]]])
    Term(GOTERM[[goids[1]]])
    Synonym(GOTERM[[goids[1]]])
    Secondary(GOTERM[[goids[1]]])
    Definition(GOTERM[[goids[1]]])
    Ontology(GOTERM[[goids[1]]])
  }
}
\keyword{datasets}

