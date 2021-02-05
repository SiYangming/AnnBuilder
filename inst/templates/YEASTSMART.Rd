\name{YEASTSMART}
\alias{YEASTSMART}
\title{Map Yeast Systematic Names to SMART IDs}
\description{
  YEASTSMART is an R environment that provides mappings between yeast
  systematic names and the associated SMART identifiers.
}
\details{
  Each yeast systematic name maps to a vector of SMART identifiers.

  Mappings were based on data provided by:

  #YEASTDOMAIN#

  Package built on #DATE#
}
\references{
  SMART website: \url{http://smart.embl-heidelberg.de/}
  and
  Data source: \url{#URL#}
}
\examples{
        # Convert to a list
        xxx <- as.list(YEASTSMART)
        # randomly display 10 probes
        sample(xxx, 10)
}
\keyword{datasets}

