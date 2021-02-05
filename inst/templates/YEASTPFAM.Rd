\name{YEASTPFAM}
\alias{YEASTPFAM}
\title{Map Yeast Systematic Names to Pfam IDs}
\description{
  YEASTPFAM is an R environment that provides mappings between yeast
  systematic names and the associated Pfam identifiers.
}
\details{
  Each yeast systematic name maps to a vector of Pfam identifiers.

  Mappings were based on data provided by:

  #YEASTDOMAIN#

  Package built on #DATE#
}
\references{
  Pfam website: \url{http://www.sanger.ac.uk/Software/Pfam/}
  and
  Data source \url{#URL#}
}
\examples{
        # Convert to a list
        xxx <- as.list(YEASTPFAM)
        # randomly display 10 probes
        sample(xxx, 10)
}
\keyword{datasets}

