\name{getUGShortName}
\alias{getUGShortName}
\alias{UGSciNames}
\alias{getShortSciName}

\title{Functions that produce short versions of organism names used by
  UniGene or for other purposes}
\description{
  From a two-word scientific name of an organism, the functions
  construct a short string used by UniGene or others to represent the
  organism.  
}
\usage{
getUGShortName(sciName)
UGSciNames()
getShortSciName(sciName)
}

\arguments{
  \item{sciName}{\code{sciName} a character string for the scientific
    name of an organism}
}
\details{
  Given a two-word scientifc name for a given organism,
  \code{\link{getUGShortName}} figures out the short version used by
  UniGene as part of the name for the file containing data for the
  organism.

  \code{\link{getShortSciName}} takes a two-word scientific name of an
  organism and returns a three-letter string begining with the first
  letter of the genus name followed by the first two letters of the
  species name.
}
\value{
  \code{\link{getUGShortName}} returns a short version of organism name
  used by UniGene.
}
\author{Jianhua Zhang}

\seealso{\code{\link{ABPkgBuilder}}}
\examples{
\dontrun{
    getUGShortName("Homo sapiens")
    getShortSciName("Homo sapiens")
}
}
\keyword{misc}


