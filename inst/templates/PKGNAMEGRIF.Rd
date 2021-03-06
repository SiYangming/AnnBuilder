\name{#PKGNAME#GRIF}
\alias{#PKGNAME#GRIF}
\title{Map between Manufacturer Identifiers and GeneRIF}
\description{
 #PKGNAME#GRIF is an R environment that provides mappings between
 manufacturer identifiers and GeneRIF identifiers.
}
\details{
  Each manufacturer identifier is mapped to a named vector containing GeneRIF
  identifiers. The name associated with each vector corresponds to the
  manufacturer identifier.

  GeneRIF is a manual method for scientists to annotate genes described
  in Entrez Gene with a phrase describing the function and a published
  paper (from PubMed) describing that function.
  
  Mappings were based on data provided by: 

  #LLSOURCE#

  Package built: #DATE#
}
\references{
\url{http://www.ncbi.nlm.nih.gov/projects/GeneRIF/GeneRIFhelp.html}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#GRIF)
    # Remove probe identifiers that do not map to any PubMed id
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # The unique PubMed identifiers for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    }
}
\keyword{datasets}

