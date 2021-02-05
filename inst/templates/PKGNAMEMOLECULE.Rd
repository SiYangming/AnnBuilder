\name{#PKGNAME##DSOURCE#MOLECULE}
\alias{#PKGNAME##DSOURCE#MOLECULE}
\title{An annotation data file that maps cMAP molecule identifiers to
  data concerning the molecule}
\description{
  #PKGNAME##DSOURCE#MOLECULE maps cMAP (NCICB Pathway Interaction Database)
  molecule identifiers to data about the molecules   
}
\details{
  This is an environment object containing key and value pairs. Keys are
  molecule identifiers and values are lists of vectors and sublists. Each
  molecule identifier is mapped to a list that has a type, extid, component, and
  member element.

  The type element is a character string describing what type of
  molecule the key molecule identifier corresponds to. Potential values for type
  include "protein", "complex", "compound", and "rna".

  The extid element is a named vector of character strings. values of
  the vector are external identifiers corresponding to the key molecule
  identifier. Names of the vector are abbreviations of the external public data
  sources from which the external identifiers are obtained. Potential values for
  vector names include "LL" (Entrez Gene identifier), "CA" (CAChemical
  Abstract number), "GO" (Gene Ontology), "KG" (KEGG), etc.

  The component element is a sublist with an identifier (molecule id
  of the component molecule), location (a GO Cellular Component
  vocabulary indicating the location of the molecule), and activity
  (an abstract term that can be "inactive", "active", "active1",
  "active2") elements. The component element only applies to complex
  molecules. Each complex molecule has a sequence of component molecules
  represented by elements of the sublist. 

  The member element is a named vector of molecule identifiers for molecules
  that belong to the same protein family as the key molecule identifier. Names
  of the vector are the key molecule identifier. 

  Mappings were based on data provided by:

  #CMAPSOURCE#

  Package built: #DATE#
}
\references{
cMAP \url{http://cmap.nci.nih.gov/PW}
}
\examples{
    require(#PKGNAME#) || stop(paste(#PKGNAME#, "unavailable"))
    xx <- as.list(#PKGNAME##DSOURCE#MOLECULE)
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

