\name{#PKGNAME##DSOURCE#PATHWAY}
\alias{#PKGNAME##DSOURCE#PATHWAY}
\title{An annotation data file that maps cMAP short pathway names to
  descriptive data about the pathway}
\description{
  #PKGNAME##DSOURCE#PATHWAY maps cMAP (NCICB Pathway Interaction Database)
  abbreviated pathway names to descriptive data about the pathway   
}
\details{
  This is an environment object containing key and value pairs. Keys are
  short pathway names and values are lists of vectors. Each
  pathway name is mapped to a list that has an identifier. name, organism,
  source, and component element.

  The identifier element is an integer for NCICB Pathway Interaction Database
  identifier.

  The name element is the full length of textual descriptive name for
  the pathway.

  The organism is a character string for an abbreviation of organism
  name (e. g. Hs = human)

  The source element is a character string indicating whether the
  interaction is a BioCarta or KEGG pathway.

  The type element is a character string describing what type of
  molecule the key molecule identifier corresponds to. Potential values for type
  include "protein", "complex", "compound", and "rna".

  The component element is a vector of integers for the identifiers of molecules
  that are involved in the pathway.  

  Mappings were based on data provided by:

  #CMAPSOURCE#

  Package built: #DATE#
}
\references{
cMAP \url{http://cmap.nci.nih.gov/PW}
}
\examples{
    require(#PKGNAME#) || stop(paste(#PKGNAME#, "unavailable"))
    xx <- as.list(#PKGNAME##DSOURCE#PATHWAY)
    if(length(xx) > 0){
        # Get the value of the first key
        xx[[1]]
        # Get the values for multiget for a few keys
        if(length(xx) >= 3){
   	    xx[1:3]
        }
    }
}
\keyword{datasets}

