\name{#PKGNAME#HGID}
\alias{#PKGNAME#HGID}
\title{A vector containing all the HomoloGeneIDs for #ORGANISM#}
\description{
  #PKGNAME#HGID is a vector containing all the HomoloGeneIDs for
  #ORGANISM#. A HomoloGeneID is a unique identifier used by
  NCBI to represent a sequence in a given organism 
}
\details{
  The HomoloGeneIDs can be used to obtain data for all the homologous
  genes in other organisms using the environment #PKGNAME#DATA.

  #PKGNAME#DATA is an R environment that maps HomoloGeneIds from
  #ORGANISM# to other species. Consult the #PKGNAME#DATA help page for
  more information.
      
  Mappings contained in homology#ORGCODE# were based on data provided
  by HomoloGene.
  
  #HGSOURCE#

  Packagebuilt #DATE#
}
\references{
\url{http://www.ncbi.nih.gov/entrez/query.fcgi?db=homologene}
}
\examples{
    require("annotate") || stop("annotate unavailable")
    get(#PKGNAME#HGID[1], #PKGNAME#DATA)
}
\keyword{datasets}

