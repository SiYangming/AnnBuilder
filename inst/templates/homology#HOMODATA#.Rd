\name{homology#HOMODATA#}
\alias{homology#HOMODATA#}
\title{Maps HomoloGeneIDs to homology data of other organisms that have
  hology relationships with genes represented by the HomoloGeneIDs
}
\description{
  Genes that share an arbitrary threshold level of similarity determined
  by alignment of matching bases are termed
  homologous. homology#HOMODATA# contains mappings between internal 
  HomologGeneIDs (used by NCBI to represent sequences of a given organism)
  and the data concerning genes in other organisms that have been
  identified to have homology reltationships with genes represented by
  the HomologGeneIDs  
}
\details{
  homology#HOMODATA# is an environment object with key and value pairs.
  Keys are HomoloGeneIDs and values are lists of homoData
  objects for other organisms that have been identified to have homology
  relationships with the genes represented by the HomoloGeneIDs. The
  list for a given internal HomologGeneID may have one to many homoData
  objects depending on whether homolgy relationships have been
  identified in only one or more other organisms.

  A homoData object has the following slots:

  homoOrg - a named vector of a single character string whose value is
  the scientific name of the organism and name the numeric code used by
  NCBI for the organism.
  
  homoLL - an integer for LocusLink id.

  homoHGID - an integer for internal HomologGeneID.

  homoACC - a character string for GenBank accession number of the best
  matching sequence of the organism.

  homoType - a single letter for the type of similarity measurement
  between the homologous genes. homoType can be either B (reciprocal
  best best between three or more organisms), b (reciprocal best match
  between two organisms), or c (curated homology relationship between
  two organisms).

  homoPS - a percentage value measured as the percent of identity of
  base pair alignment between the homologous sequences. 

  homoURL - a url to the source if the homology relationship is a curated
  orthology.

  HomoData objects with homoType = B or b will not have any
  value for homoURL and objects with homoType = c will not have any
  value for homoPS. 

  The value of a given slot of a homoData object can be accessed by
  using homoOrg(xxx) for the homoOrg slot of the homoData object xxx,
  for example.  

  Mappings contained in homology#HOMODATA# were based on data provided
  by HomoloGene.
  
  #HGSOURCE#

  Packagebuilt #DATE#
}
\references{
\url{ftp://ftp.ncbi.nih.gov/pub/HomoloGene/hmlg.ftp}
}
\examples{
    require("annotate") || stop("annotate unavailable")
    # Convert homology#HOMODATA# to a list
    xx <- as.list(homology#HOMODATA#)

    if(length(xx) > 0){
        # Gets the homology data for the first key of homology#HOMODATA#
        homodata <- xx[[1]]
        # Gets the organism name, LocusLink id, internal HomoloGene id,
        # GenBank accession number, type of similarity, percent
        # similarity, and URL for percent similarity for the first
        # homoData object in the list obtained previously
        homoOrg(homodata[[1]])
        homoLL(homodata[[1]])
        homoHGID(homodata[[1]])
        homoACC(homodata[[1]])
        homoType(homodata[[1]])
        homoPS(homodata[[1]])
        homoURL(homodata[[1]])
        # Gets the values for more than one keys
        xx[1:3]
    }
}
\keyword{datasets}

