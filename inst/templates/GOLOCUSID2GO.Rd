\name{GOLOCUSID2GO}
\alias{GOLOCUSID2GO}
\title{Entrez Gene to Gene Ontology (GO) mapping}
\description{
  This data set describes associations between Entrez Gene identifiers
  and GO identifiers as reported by NCBI.  The format is an R
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
  supports the association of the GO id to the Entrez Gene id. The
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

  Entrez Gene ids for which no GO associations exist are left out of the
  environment.

  Mappings were based on data provided by:

  #LLSOURCE#

  Package built: #DATE#
}
\references{
  \url{http://www.ncbi.nih.gov/entrez/query.fcgi?db=gene}
}
\examples{
library("GO")
## Select an Entrez Gene id
egId <- "10043"
egIdGoData <- GOLOCUSID2GO[[egId]]

## Display GO ids associated with the Entrez Gene id
sapply(egIdGoData, function(x) x$GOID)

## The ontology categories of the GO ids
sapply(egIdGoData, function(x) x$Ontology)

## The evidence codes for the GO ids
sapply(egIdGoData, function(x) x$Evidence)
}
\keyword{datasets}

