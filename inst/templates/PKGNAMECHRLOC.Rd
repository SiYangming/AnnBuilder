\name{#PKGNAME#CHRLOC}
\alias{#PKGNAME#CHRLOC}
\title{Map Manufacturer IDs to Chromosomal Location}
\description{
  #PKGNAME#CHRLOC is an R environment that maps manufacturer identifiers to the
  starting position of the gene. The position of a gene is
  measured as the number of base pairs.
}
\details{
  Each manufacturer identifier maps to a named vector of chromosomal locations,
  where the name indicates the chromosome. Due to inconsistencies that
  may exist at the time the environment was built, these vectors may
  contain more than one chromosome and/or location. If the chromosomal
  location is unknown, the vector will contain an \code{NA}.
  
  Chromosomal locations on both the sense and antisense strands are
  measured as the number of base pairs from the p (5' end of the sense
  strand) to q (3' end of the sense strand) arms. Chromosomal locations on the
  antisense strand have a leading "-" sign (e. g. -1234567). 
  
  Package built #DATE#
}
\examples{
    # Covert to a list
    xx <- as.list(#PKGNAME#CHRLOC)
    # Remove probe identifiers that do not map to any CHRLOC
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # Gets the location for the first five probes
        xx[1:5]
        # Gets the first one
        xx[[1]]
    }
}
\keyword{datasets}

