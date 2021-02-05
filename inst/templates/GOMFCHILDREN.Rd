\name{GOMFCHILDREN}
\alias{GOMFCHILDREN}
\title{Annotation of GO Identifiers to their Molecular Function Children} 
\description{
  This data set describes associations between GO molecular function (MF)
  terms and their direct children MF terms, based on the directed acyclic
  graph (DAG) defined by the Gene Ontology Consortium. The format is an R
  environment mapping the GO MF terms to all direct children terms, where a
  direct child term is a more specific GO term that is immediately preceded
  by the given GO term in the DAG.
}
\details{
  Each GO MF term is mapped to a vector of children GO MF terms.

  Molecular function is defined as the tasks performed by individual
  gene products; examples are transcription factor and DNA helicase as
  defined by Gene Ontology Consortium.  
  
  Mappings were based on data provided by:

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
    xx <- as.list(GOMFCHILDREN)
    # Remove GO ids that do not have any children
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # Get the children GO ids for the first elents of xx
        goids <- xx[[1]]
        # Find out the GO terms for the first parent goid
        GOID(GOTERM[[goids[1]]])
        Term(GOTERM[[goids[1]]])
        Synonym(GOTERM[[goids[1]]])
        Secondary(GOTERM[[goids[1]]])
        Definition(GOTERM[[goids[1]]])
        Category(GOTERM[[goids[1]]])
    }
}
\keyword{datasets}

