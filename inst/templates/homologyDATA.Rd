\name{#PKGNAME#DATA}
\alias{#PKGNAME#DATA}
\title{Genes that are Homologous to #ORGANISM# 
}
\description{
  This data set gives mappings between #ORGANISM# genes and other
  species based on homology. Genes that share an arbitrary threshold
  level of similarity determined by alignment of matching bases are termed
  homologous. The format is an R environment that maps HomoloGeneIds
  (which are used by NCBI to represent sequences of an organism) to a
  list of organisms that have homologous sequences.
}
\details{
  Each HomoloGeneId is mapped to a list of lists. The names of the outer
  list are the numeric code used for that organism by the NCBI. This is
  included for convenience when processing the list using \code{lapply}.

  Each inner list consists of five named elements; homoORG, homoType,
  homoPS, homoHGID, and homoURL.

  The homoORG element indicates the organism for which there is a
  homologous gene. This will be an abbreviation of the first letter of
  the genus and the first two letters of the species. For instance, Homo
  sapiens is abbreviated \emph{hsa}.

  The homoType element indicates the type of similarity
  measurement. This can be either \emph{B}, indicating a reciprocal best
  best match between three or more organisms, \emph{b}, indicating a
  reciprocal best match between two organisms, or \emph{c}, indicating a
  curated homology relationship between two organisms.

  The homoPS element gives the percentage of identical base pair
  alignments between the homologous genes.

  The homoHGID element gives the internal HomoloGeneId.
  
  The homoURL element gives a url to the source for a curated homology.

  Note that there will not be a value for homoURL if the homoType is B
  or b, and there will not be a value for homoPS if the homoType is c.

  Mappings contained in homologyDATA were based on data provided
  by HomoloGene.
  
  #HGSOURCE#

  Packagebuilt #DATE#
}
\references{
\url{ftp://ftp.ncbi.nih.gov/pub/HomoloGene/hmlg.ftp}
}
\examples{
    require("annotate") || stop("annotate unavailable")
    # Convert homologyDATA to a list
    xx <- as.list(#PKGNAME#DATA)

    if(length(xx) > 0){
        # Gets the homology data for the first key of #PKGNAME#DATA
        homodata <- xx[[1]]
        # Gets the organism name, LocusLink id, internal HomoloGene id,
        # GenBank accession number, type of similarity, percent
        # similarity, and URL for percent similarity for the first
        # homoData object in the list obtained previously
        homodata[[1]][["homoOrg"]]
        homodata[[1]][["homoHGID"]]
        homodata[[1]][["homoType"]]
        homodata[[1]][["homoPS"]]
        homodata[[1]][["homoURL"]]
        # Gets the values for more than one keys
        xx[1:3]
    }
}
\keyword{datasets}

