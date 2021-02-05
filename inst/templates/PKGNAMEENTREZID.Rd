\name{#PKGNAME#ENTREZID}
\alias{#PKGNAME#ENTREZID}
\alias{#PKGNAME#LOCUSID}
\title{Map between Manufacturer Identifiers and Entrez Gene}
\description{
  #PKGNAME#ENTREZID is an R environment that provides mappings between
  manufacturer identifiers and Entrez Gene identifiers.
}
\details{
  Each manufacturer identifier is mapped to a vector of Entrez Gene identifiers. An
  \code{NA} is assigned to those manufacturer identifiers that can not be mapped
  to an Entrez Gene identifier at this time.
  
  If a given manufacturer identifier can be mapped to different Entrez Gene identifiers
  from various sources, we attempt to select the common identifiers. If a
  concensus cannot be determined, we select the smallest identifier.
  
  Mappings were based on data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#ENTREZID)
    # Remove probe identifiers that do not map to any ENTREZID
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # The ENTREZIDs for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    }
}

\keyword{datasets}

