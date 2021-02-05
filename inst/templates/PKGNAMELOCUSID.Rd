\name{#PKGNAME#LOCUSID}
\alias{#PKGNAME#LOCUSID}
\title{Map between Manufacturer Identifiers and Entrez Gene}
\description{
  #PKGNAME#LOCUSID is an R environment that provides mappings between
  manufacturer ids and Entrez Gene ids.
}
\details{
  Each manufacturer id is mapped to a vector of Entrez Gene ids. An
  \code{NA} is assigned to those manufacturer ids that can not be mapped
  to an Entrez Gene id at this time.
  
  If a given manufacturer id can be mapped to different Entrez Gene ids
  from various sources, we attempt to select the common ids. If a
  concensus cannot be determined, we select the smallest id.
  
  Mappings were based on data provided by:

  #LLSOURCE#
  
  Package built: #DATE#
}
\references{
\url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#LOCUSID)
    # Remove probe ids that do not map to any LOCUSID
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # The LOCUSIDs for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    }
}

\keyword{datasets}

