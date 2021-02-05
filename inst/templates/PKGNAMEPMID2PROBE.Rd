\name{#PKGNAME#PMID2PROBE}
\alias{#PKGNAME#PMID2PROBE}
\title{Map between PubMed Identifiers and Manufacturer Identifiers}
\description{
  #PKGNAME#PMID2PROBE is an R environment that provides mappings between
  PubMed ids and manufacturer ids.
}
\details{
  Each PubMed id is mapped to a named vector of manufacturer ids. The
  name represents the PubMed id and the vector contains all manufacturer
  ids that are represented by that PubMed id. The length of the vector
  may be one or longer, depending on how many manufacturer ids are
  mapped to a given PubMed id.

  Titles, abstracts, and possibly full texts of articles can be obtained
  from PubMed by providing a valid PubMed id. The pubmed function of
  annotate can also be used for the same purpose
  
  Mappings were based on data provided by:
  
  #LLSOURCE# 
  
  Package built #DATE#
}

\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#PMID2PROBE)
    if(length(xx) > 0){
        # The probe ids for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
        if(interactive() && require(annotate)){
            # Gets article information as XML files for a PubMed id
            xmls <- pubmed(names(xx)[1], disp = "data")
            # Views article information using a browser
            pubmed(names(xx)[1], disp = "browser")
        }
    }
}

\keyword{datasets}

