\name{#PKGNAME#PATH}
\alias{#PKGNAME#PATH}
\title{Mappings between probe identifiers and KEGG pathway identifiers}
\description{
  KEGG (Kyoto Encyclopedia of Genes and Genomes) maintains pathway data
  for various organisms. #PKGNAME#PATH maps probe identifiers to the
  identifiers used by KEGG for pathways in which the genes represented
  by the probe identifiers are involved
}
\details{
  Each KEGG pathway has a name and identifier. Pathway name for a given pathway
  identifier can be obtained using the KEGG data package that can either be
  built using AnnBuilder or downloaded from Bioconductor
  \url{http://www.bioconductor.org}.

  Graphic presentations of pathways are searchable at
  url{http://www.genome.ad.jp/kegg/pathway.html} by using pathway identifiers as
  keys.   

  Mappings were based on data provided by:

  #KEGGSOURCE#

  #LLSOURCE#

  Package built #DATE#
}
\references{
  \url{http://www.genome.ad.jp/kegg/}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#PATH)
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

