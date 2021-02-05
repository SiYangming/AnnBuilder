\name{addNamespace}
\alias{addNamespace}
\alias{sealEnvs}
\title{Functions to add namespaces for data files or seal the
  environment objects in the data subdirectory} 
\description{
  Given the name of a data package and the path, the functions add
  namespaces for data files in the data subdirectory or seal the
  environment objects in the data subdirectory
}
\usage{
addNamespace(pkgName, pkgPath, hidePattern = c("QC", "MAPCOUNTS"))
sealEnvs(pkgName, pkgPath)
}

\arguments{
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package whose data subdirectory contains data files to be put
    in the NAMESPACE}
  \item{pkgPath}{\code{pkgPath} a character stirng for the path where a
    data package of interest resides}
  \item{hidePattern}{\code{hidePattern} a vector of character strings
    whose patterns match the data files in the data subdirectory that
    will not included in the NAMESPACE}
}
\details{
  These functions are mainly for manipulating data files for annotation
  data packages and may not be of other usages.
}
\value{
  The functions returns invisible()
}
\author{Jianhua Zhang}

\examples{
  # No examples provided
}
\keyword{misc}
