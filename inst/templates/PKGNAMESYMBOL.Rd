\name{#PKGNAME#SYMBOL}
\alias{#PKGNAME#SYMBOL}
\title{Map between Manufacturer Identifiers and Gene Symbols}
\description{
  #PKGNAME#SYMBOL is an R environment that provides mappings between
  manufacturer ids and gene abbreviations. 
}
\details{
  Each manufacturer id is mapped to an abbreviation for the
  corresponding gene. An \code{NA} is reported if there is no known
  abbreviation for a given gene.
  
  Symbols typically consist of 3 letters that define either a single
  gene (ABC) or multiple genes (ABC1, ABC2, ABC3). Gene symbols can be
  used as key words to query public databases such as Entrez Gene.
  
  Mappings were based on data provided by:
  
  #LLSOURCE#
  
  Package built #DATE#
}
\references{
  \url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#SYMBOL)
    if(length(xx) > 0){
        # The symbols for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    }
}

\keyword{datasets}

