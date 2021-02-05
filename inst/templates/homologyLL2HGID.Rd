\name{#PKGNAME#LL2HGID}
\alias{#PKGNAME#LL2HGID}
\alias{#PKGNAME#ACC2HGID}
\alias{#PKGNAME#HGID2LL}
\alias{#PKGNAME#HGID2ACC} 
\title{Mappings between public database identifiers and NCBI's Internal
  HomoloGeneIDs for #ORGANISM#}  
\description{
  R environments that provides mappings between public database identifiers and internal
  HomoloGeneIDs (HGID) used by NCBI to represent sequences for
  #ORGANISM#. Public database identifiers include Entrez Gene identifiers (LL) and GenBank
  Accesion/RefSeq numbers (ACC).     
}
\details{
  Each R environment provides mapping between either HomoloGeneIds and
  database identifiers. or the converse. The mapping direction is indicated by
  the name of the environment. For instance, #PKGNAME#ACC2HGID provides
  mappings between GenBank accession numbers and HomoloGeneIds.

  For each environment, a given identifier will map to a named vector of identifiers.
  the names of which are the NCBI code for that organism.
 
  Mappings contained were based on data provided by HomoloGene.
  
  #HGSOURCE#

  Packagebuilt #DATE#
}
\references{
\url{http://www.ncbi.nih.gov/entrez/query.fcgi?db=homologene}
}
\examples{
    require("annotate") || stop("annotate unavailable")
    xx <- as.list(#PKGNAME#LL2HGID)
        if(length(xx) > 5){
                # Get the value of the first key
                xx[2]
                # Get the value for a few keys
                if(length(xx) >= 3){
                        xx[2:4]
                }
        }
}
\keyword{datasets} 

