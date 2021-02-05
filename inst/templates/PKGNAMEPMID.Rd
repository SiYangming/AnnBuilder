\name{#PKGNAME#PMID}
\alias{#PKGNAME#PMID}
\title{Map between Manufacturer Identifiers and PubMed Identifiers}
\description{
   #PKGNAME#PMID is an R environment that provides mappings between
   manufacturer identifiers and PubMed identifiers.
}
\details{
  Each manufacturer identifier is mapped to a named vector of PubMed identifiers. The
  name associated with each vector corresponds to the manufacturer
  identifier. The length of the vector may be one or greater, depending on how
  many PubMed identifiers a given manufacturer identifier is mapped to. An \code{NA} is
  reported for any manufacturer identifier that cannot be mapped to a PubMed identifier.

  Titles, abstracts, and possibly full texts of articles can be obtained
  from PubMed by providing a valid PubMed identifier. The pubmed function of
  annotate can also be used for the same purpose.

  Mappings were based on data provided by:

  #LLSOURCE#

  Package built #DATE#
}
\references{
  \url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#PMID)
    if(length(xx) > 0){
        # The probe identifiers for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
        if(interactive() && !is.null(xx[[1]]) && !is.na(xx[[1]])
           && require(annotate)){
            # Gets article information as XML files
            xmls <- pubmed(xx[[1]], disp = "data")
            # Views article information using a browser
            pubmed(xx[[1]], disp = "browser")
        }
    }
}
\keyword{datasets}

