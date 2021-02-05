\name{#PKGNAME#GENENAME}
\alias{#PKGNAME#GENENAME}
\title{Map Between Manufacturer IDs and Genes}
\description{
  #PKGNAME#GENENAME is an R environment that maps manufacturer
  identifiers to the corresponding gene name.
}
\details{
  Each manufacturer id maps to a named vector containing the gene
  name. The vector name corresponds to the manufacturer id. If the gene
  name is unknown, the vector will contain an \code{NA}.
  
  Gene names currently include both the official (validated by a
  nomenclature committee) and preferred names (interim selected for
  display) for genes. Efforts are being made to differentiate the two
  by adding a name to the vector.

  Mappings were based on data provided by:

  #LLSOURCE#
  
  Package built #DATE#
}

\examples{
     # Convert to a list
     xx <- as.list(#PKGNAME#GENENAME)
     # Remove probes that do not map to any GENENAME
     xx <- xx[!is.na(xx)]
     if(length(xx) > 0){
         # Gets the gene names for the first five probe ids
         xx[1:5]
         # Get the first one
         xx[[1]]
     }
}
\keyword{datasets}

