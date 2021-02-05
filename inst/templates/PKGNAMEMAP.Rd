\name{#PKGNAME#MAP}
\alias{#PKGNAME#MAP}
\title{Map between Manufacturer Identifiers and cytogenetic maps/bands}
\description{
   #PKGNAME#MAP is an R environment that provides mappings between
   manufacturer ids and cytoband locations.
}
\details{
  Each manufacturer id is mapped to a vector of cytoband locations. The
  vector length may be one or longer, if there are multiple reported
  chromosomal locations for a given gene. An \code{NA} is reported for
  any manufacturer ids that cannot be mapped to a cytoband at this time.
  
  Cytogenetic bands for most higher organisms are labeled p1, p2, p3,
  q1, q2, q3 (p and q are the p and q arms), etc., counting from the
  centromere out toward the telomeres. At higher resolutions, sub-bands
  can be seen within the bands. The sub-bands are also numbered from the
  centromere out toward the telomere. Thus, a label of 7q31.2 indicates
  that the band is on chromosome 7, q arm, band 3, sub-band 1, and
  sub-sub-band 2.   

  The physical location of each band on a chromosome can be obtained
  from another environment named "organism"CYTOLOC in a separate data
  package for human(humanCHRLOC), mouse(mouseCHRLOC), and rat(ratCHRLOC).  
  
  Mappings were based on data provided by:

  #LLSOURCE#

  Package built #DATE#
}
\references{
  \url{http://www.ncbi.nlm.nih.gov}
}
\examples{
    # Convert the environment to a list
    xx <- as.list(#PKGNAME#MAP)
    # Remove probe ids that do not map to any cytoband
    xx <- xx[!is.na(xx)]
    if(length(xx) > 0){
        # The cytobands for the first two elements of XX
        xx[1:2]
        # Get the first one
        xx[[1]]
    }
}
\keyword{datasets}

