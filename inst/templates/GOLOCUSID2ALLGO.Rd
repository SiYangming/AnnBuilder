\name{GOLOCUSID2ALLGO}
\alias{GOLOCUSID2ALLGO}
\title{Entrez Gene to Gene Ontology (GO) Mapping with Child Nodes} 
\description{
  This data set describes associations between Entrez Gene identifiers
  and GO identifiers as reported by NCBI as well as associations derived
  from the parent-child relations from GO.  The format is an R
  environment mapping Entrez Gene ids to a list of GO ids that are
  associated with it.  In addition to the GO id, the GO ontology
  category and the evidence code (the reason for the association) is
  given.  
}
\details{
  Each Entrez Gene id is mapped to a list of lists.  The names on the
  outer list are GO ids.  Each inner list consists of three named
  elements: GOID, Ontology, and Evidence.

  The GOID element matches the GO id named in the outer list and is
  included for convenience when processing the data using \code{lapply}.

  The Ontology element indicates which of the three Gene Ontology
  categories this identifier belongs to.  The categories are biological
  process (BP), cellular component (CC), and molecular function (MF).

  The Evidence element contains a code indicating what kind of evidence
  supports the association of the GO id to the LocusLink id. The
  evidence codes in use include:


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
  
  NA is assigned to Entrez Gene id that can not be mapped to any GO id at
  this time.

  Mappings were based on data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
    \url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} 
}
\examples{
  	require("GO") || stop("GO unavailable")
        # Convert the environment object to a list
	xx <- as.list(GOLOCUSID2ALLGO)
        # Remove Entrez Gene ids that are not mapped to any GO id
        xx <- xx[!is.na(xx)]
	if(length(xx) > 0){
		# Get the GO ids for the first two elents of xx
		goids  <- xx[1:2]
                # Get the evidence code goids
                evds <- sapply(goids, names)
	}
}
\keyword{datasets}




