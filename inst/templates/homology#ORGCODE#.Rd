\name{homology#ORGCODE#}
\alias{homology#ORGCODE#}
\title{A vector containing all the HomoloGeneIDs for #ORGANISM#}
\description{
  homology#ORGCODE# is a vector containing all the HomoloGeneIDs for
  #ORGANISM#. A HomoloGeneID is a unique identifier used by
  NCBI to represent a sequence in a given organism 
}
\details{
  The HomoloGeneIDs can be used to obtain data for all the homologous
  genes in other organisms that have been identified to have homology
  relationships with the sequences represented by the HomoloGeneIDs
  using the environment homologyDATA.

  homologyDATA is an environment object with key and value pairs.
  Keys are HomoloGeneIDs and values are lists of homoData
  objects for other organisms that have been identified to have homology
  relationships with the genes represented by the HomoloGeneIDs.
      
  Mappings contained in homology#ORGCODE# were based on data provided
  by HomoloGene.
  
  #HGSOURCE#

  Packagebuilt #DATE#
}
\references{
\url{http://www.ncbi.nih.gov/entrez/query.fcgi?db=homologene}
}
\examples{
    require("annotate") || stop("annotate unavailable")
    get(homology#ORGCODE#[1], homologyDATA)
}
\keyword{datasets}

