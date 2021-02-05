\name{#PKGNAME#ORGCODE}
\alias{#PKGNAME#ORGCODE}
\title{A data frame for mappings between organism scientific names, NCBI
  code, and three-letter short names for various organisms}
\description{
  This is a data frame that maps the scientific names of various
  organisms to the numeric code used by NCBI and three-letter short
  names of the organisms.
}
\details{ 
  Data provided by HomoloGene.
  
  #HGSOURCE#

  Packagebuilt #DATE#
}
\references{
\url{http://www.ncbi.nih.gov/entrez/query.fcgi?db=homologene}
}
\examples{
    require("annotate") || stop("annotate unavailable")
    #PKGNAME#ORGCODE
}
\keyword{datasets}

