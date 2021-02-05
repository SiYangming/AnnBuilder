\name{#PKGNAME#PATH2PROBE}
\alias{#PKGNAME#PATH2PROBE}
\title{Map between Kyoto Encyclopedia of Genes and Genomes (KEGG)
  pathway identifiers and Manufacturer Identifiers}
\description{
 #PKGNAME#PATH2PROBE is an R environment that provides mappings between
 KEGG ids and manufacturer ids.
}
\details{
  Each KEGG id is mapped to a named vector of manufacturer ids. The name
  represents the KEGG id and the vector contains all manufacturer ids
  that are found in that particular pathway. An \code{NA} is reported
  for any KEGG id that cannot be mapped to any manufacturer ids.

  Pathway name for a given pathway id can be obtained using the KEGG
  data package that can either be built using AnnBuilder or downloaded
  from Bioconductor \url{http://www.bioconductor.org}.

  Graphic presentations of pathways are searchable at
  \url{http://www.genome.ad.jp/kegg/pathway.html} using pathway ids as keys.
  
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
    xx <- as.list(#PKGNAME#PATH2PROBE)
    # Remove pathway ids that do not map to any probe id
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # The probe ids for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    } 
}

\keyword{datasets}

