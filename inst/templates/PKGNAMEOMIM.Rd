\name{#PKGNAME#OMIM}
\alias{#PKGNAME#OMIM}
\title{Map between Manufacturer Identifiers and Mendelian Inheritance in
  Man (MIM) identifiers}
\description{
  #PKGNAME#OMIM is an R environment that provides mappings between
  manufacturer ids and OMIM ids.
}
\details{
  Each manufacturer id is mapped to a vector of OMIM ids. The vector
  length may be one or longer, depending on how many OMIM ids the
  manufacturer id maps to. An \code{NA} is reported for any manufacturer
  id that cannot be mapped to an OMIM id at this time.
  
  OMIM is based upon the book Mendelian Inheritance in Man
  (V. A. McKusick) and focuses primarily on inherited or heritable
  genetic diseases. It contains textual information, pictures, and reference
  information that can be searched using various terms, among which the
  MIM number is one. 
  
  Mappings were based on data provided by:

  #LLSOURCE#

  Package built #DATE#
}
\references{
  \url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene}
  \url{http://www3.ncbi.nlm.nih.gov/entrez/query.fcgi?db=OMIM}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#OMIM)
    # Remove probe ids that do not map to any MIM number
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # The MIM numbers for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    }
}

\keyword{datasets}

