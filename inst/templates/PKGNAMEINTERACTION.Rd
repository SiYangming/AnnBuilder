\name{#PKGNAME##DSOURCE#INTERACTION}
\alias{#PKGNAME##DSOURCE#INTERACTION}
\title{An annotation data file that maps cMAP molecule identifiers to
  data concerning the interactions between molecules}
\description{
  #PKGNAME##DSOURCE#INTERACTION maps cMAP (NCICB Pathway Interaction Database)
  molecule identifiers to data about the interactions of the molecule
  represented by the identifiers and other molecules   
}
\details{
  This is an environment object containing key and value pairs. Keys are
  molecule identifiers and values are lists of vectors and sublists. Each
  molecule identifier is mapped to a list that has a source, process,
  reversible, condition, and component element.

  The source element is a character string indicating whether the
  interactions between molecules are from a BioCarta or KEGG
  pathway. 

  The process is a character string describing the process the key
  molecule is involved. Potential values include "reaction",
  "modification", "transcription", "translocation", "macroprocess" or a
  more specific subtype of macroprocess including any term from the GO
  Biological Process vocabulary. 

  The reversible element is a boolean indicating whether the interaction
  is reversible. 

  The condition element is a character string indicating the biological
  process the interactions take place. Potential values include any term
  from the GO Biological Process vocabulary.
  
  The component element contains sublists of vectors. Each key molecule
  identifier has a sequence of component represented by elements of the
  sublist. Each sublist has an identifier (molecule identifier of the interacting
  molecule), edge (indicating the way two molecule interact. Potential
  values include "input", "agent", "inhibitor", and "output"), role (the
  function of the key molecule. Potential values include any term from
  the GO Molecular Function vocabulary), location (a GO Cellular Component
  vocabulary indicating the location of the interaction), and activity
  (an abstract term that can be "inactive", "active", "active1",
  "active2") elements.   

  Mappings were based on data provided by:

  #CMAPSOURCE#

  Package built: #DATE#
}
\references{
cMAP \url{http://cmap.nci.nih.gov/PW}
}
\examples{
    require(#PKGNAME#) || stop(paste(#PKGNAME#, "unavailable"))
    xx <- as.list(#PKGNAME##DSOURCE#INTERACTION)
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

