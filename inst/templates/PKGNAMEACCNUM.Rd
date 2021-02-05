\name{#PKGNAME#ACCNUM}
\alias{#PKGNAME#ACCNUM}
\title{Map Manufacturer ids to GenBank Accession Numbers}  
\description{ 
  #PKGNAME#ACCNUM is an R environment that contains mappings between a
  manufacturer's ids and GenBank accession numbers.    
}
\details{
  This environment was produced by first mapping the manufacturer's ids to
  Entrez Gene identifiers
  \url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene}. The
  mapped Entrez Gene ids then serve as the point of linkage to the
  GenBank accession numbers.

  Each manufacturer id maps to a vector containing a GenBank accession number.

  Package built: #DATE#
}

\examples{
        # Convert to a list
        xx <- as.list(#PKGNAME#ACCNUM)
        # Remove probe ids that do not map to any ACCNUM
        xx <- xx[!is.na(xx)]
        if(length(xx) > 0){
        # Get the ACCNUM for the first five probes
        xx[1:5]
        # Get the first one
        xx[[1]]
}
}
\keyword{datasets}

