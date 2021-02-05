\name{#PKGNAME#SUMFUNC}
\alias{#PKGNAME#SUMFUNC}
\title{Map between Manufacturer Identifiers and Gene Function Summaries}
\description{
  #PKGNAME#SUMFUNC is an R environment that provides mappings between
  manufacturer identifiers and brief summaries of functions of the gene products.  
}
\details{
  Each manufacturer identifier is mapped to a vector containing a brief summary
  of the function for the gene product. This information was previously
  supplied by Entrez Gene, but with the change to Entrez Gene these data
  may no longer be available, and the environment may be discontinued.
  
  Mappings were based on data provided by:

  #LLSOURCE#

  Package built #DATE#
}
\references{
   \url{http://www.ncbi.nlm.nih.gov/LocusLink}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#SUMFUNC)
    if(length(xx) > 0){
        # The function summaries for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    }
}

\keyword{datasets}

