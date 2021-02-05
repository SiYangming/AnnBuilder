\name{GOLOCUSID}
\alias{GOLOCUSID}
\title{Gene Ontology (GO) to Entrez Gene Mapping} 
\description{
  This data set gives mappings between GO identifiers and Entrez Gene
  identifiers as reported by NCBI. The format is an R environment mapping
  GO ids to a vector of Entrez Gene ids that are associated with it. In
  addition to the Entrez Gene id, the evidence code (the reason for the
  association) is also given.
}
\details{
  Each GO id is mapped to a named vector of Entrez Gene ids. The name
  associated with each Entrez Gene id corresponds to the evidence code
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

  NA is assigned to GO ids that can not be mapped to any Entrez Gene id at
  this time.  

  Mappings were based on data provided by:

  #GOSOURCE#

  #LLSOURCE#

  Package built: #DATE#
}
\references{
  \url{http://www.geneontology.org/} and
\url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} 
}
\examples{
    require("GO") || stop("GO unavailable")
    # Convert the environment object to a list
    xx <- as.list(GOLOCUSID)
    # Remove GO ids that are not mapped to any Entrez Gene id
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # Get the Entrez Gene ids for the first two elents of xx
        llids  <- xx[1:2]
        # Get the evidence code llids
        evds <- sapply(llids, names)
    }
}
\keyword{datasets}

