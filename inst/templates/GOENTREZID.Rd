\name{GOENTREZID}
\alias{GOENTREZID}
\alias{GOLOCUSID}
\title{Gene Ontology (GO) to Entrez Gene Mapping} 
\description{
  This data set gives mappings between GO identifiers and Entrez Gene
  identifiers as reported by NCBI. The format is an R environment mapping
  GO identifiers to a vector of Entrez Gene identifiers that are associated with it. In
  addition to the Entrez Gene identifier. the evidence code (the reason for the
  association) is also given.
}
\details{
  Each GO identifier is mapped to a named vector of Entrez Gene identifiers. The name
  associated with each Entrez Gene identifier corresponds to the evidence code
  for that GO identifier. The evidence code indicates what kind of evidence
  supports the association between the GO and Entrez Gene identifiers. Evidence
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

  NA is assigned to GO identifiers that can not be mapped to any Entrez Gene identifier at
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
    xx <- as.list(GOENTREZID)
    # Remove GO identifiers that are not mapped to any Entrez Gene id
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # Get the Entrez Gene identifiers for the first two elents of xx
        llids  <- xx[1:2]
        # Get the evidence code llids
        evds <- sapply(llids, names)
    }
}
\keyword{datasets}

