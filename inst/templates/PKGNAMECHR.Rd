\name{#PKGNAME#CHR}
\alias{#PKGNAME#CHR}
\title{Map Manufacturer IDs to Chromosomes}
\description{
 #PKGNAME#CHR is an R environment that provides mappings between a
 manufacturer identifier and the chromosome that contains the gene of interest. 
}
\details{
  Each manufacturer identifier maps to a vector of chromosomes. Due to
  inconsistencies that may exist at the time the environment was built,
  the vector may contain more than one chromosome (e.g., the identifier may map
  to more than one chromosome). If the chromosomal location is unknown,
  the vector will contain an \code{NA}.
  
  Mappings were based on data provided by:
  
  #LLSOURCE#
  
  Package built: #DATE#
}

\examples{
        # Convert to a list
        xx <- as.list(#PKGNAME#CHR)
        # Remove probe identifiers that do not map to any CHRLOC
        xx <- xx[!is.na(xx)]
        if(length(xx) > 0){
            # Get the chromosome number for the first five probes
            xx[1:5]
            # Get the first one
            xx[[1]]
        }
}
\keyword{datasets}

