\name{YEASTINTERPRO}
\alias{YEASTINTERPRO}
\title{Map Yeast Systematic Names to InterPro IDs}
\description{
  YEASTINTERPRO is an R environment that provides mappings between yeast
  systematic names and the associated InterPro identifiers.
}
\details{
  Each yeast systematic name maps to a vector of InterPro identifiers.

  Mappings were based on data provided by:

  #YEASTDOMAIN#

  Package built on #DATE#
}
\references{
  InterPro website: \url{http://www.ebi.ac.uk/interpro/}
  and
  Data source: \url{#URL#}
}
\examples{
        # Convert to a list
        xxx <- as.list(YEASTINTERPRO)
        # randomly display 10 probes
        sample(xxx, 10)
}
\keyword{datasets}

