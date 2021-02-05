\name{GOALLLOCUSID}
\alias{GOALLLOCUSID}
\title{Annotation of GO identifiers by all the corresponding Entrez Gene
  identifiers} 
\description{
This is an R environment (hash table) mapping GO ids to all the Entrez Gene
ids corresponding to the GO id and its direct or indirect
offspring. The direct or indirect offspring of a given GO id are the
nodes that are directly linked to the node for the GO id or indirectly
linked to the GO id through nodes that are linked to the GO id based on
the directed acyclic graph defined by Gene Ontology Consortium
}
\details{
  GO ids are Keys and the corresponding Entrez Gene ids are Values. Values
  are named vectors of length 1 or greater depending on whether a given
  GO id can be mapped to one or more Entrez Gene ids. Names for
  values are the evidence codes for the GO id (if an evidence code was
  provided by source data). The evidence codes that are in use include:
  
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
	xx <- as.list(GOALLLOCUSID)
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

