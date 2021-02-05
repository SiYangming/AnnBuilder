\name{#PKGNAME#ARACYC}
\alias{#PKGNAME#ARACYC}
\title{Mappings between probe identifiers and KEGG pathway identifiers}
\description{
  AraCyc \url{http://www.arabidopsis.org/tools/aracyc/} maintains pathway data
  for Arabidopsis thaliana. #PKGNAME#ARACYC maps probe identifiers to the
  common names of the pathways in which the genes represented
  by the probe identifiers are involved. Information is obtained from
  AraCyc.
}
\details{
  Package built #DATE#
}
\references{
  \url{http://www.genome.ad.jp/kegg/}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#ARACYC)
    # Remove probe identifiers that do not map to any pathway id
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # The pathway identifiers for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    } 
}
\keyword{datasets}

