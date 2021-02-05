\name{#PKGNAME#CYTOLOC}
\alias{#PKGNAME#CYTOLOC}
\title{An annotation data file for Cytoband locations on chromosomes}
\description{
  #PKGNAME#CYTOLOC maps chromosome numbers and the locations of
  cytobands on chromosoms    
}
\details{
  This is an environment object containing key and value pairs. Keys are
  chromosome numbers and values are the locations of cytobands on
  correponding chromosoms. The mapped values are lists of named
  vectors. The names of lists are cytoband ids (e. g. qA2, ...). Each
  list contains a vector of two elements of integers for the starting
  and ending locations of the band on the chromosome defined by the key
  the lists mapped to. Names of the vectors indicate whether the value
  is for the starting or ending location.
  
  Mappings were based on the following source(s):

  #GPSOURCE#

  Package built: #DATE#
}

\examples{
    bands <- as.list(#PKGNAME#CYTOLOC)
    # cytobands on chromosome number 1
    names(bands[[1]])
    # The start and end locations for one band on chromosome number 1
    bands[[1]][[1]]
}
\keyword{datasets}

