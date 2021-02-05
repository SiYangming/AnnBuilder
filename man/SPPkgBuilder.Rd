\name{SPPkgBuilder}
\alias{SPPkgBuilder}
\alias{key}
\alias{getDetailV}
\alias{getEnvNames}
\alias{isOneToOne}
\title{A function to build a data pckage using Swiss-Prot protein data}
\description{
  Given the URL to Swiss-Prot protein data, this function creates a data
  package with the data stored as R environment objects in the data
  directory 
}
\usage{
SPPkgBuilder(pkgPath, version, author, fromWeb = TRUE, url =
"ftp://ftp.ebi.ac.uk/pub/databases/swissprot/release/sprot41.dat")
getDetailV(key)
getEnvNames()
isOneToOne(envName)
}
\arguments{
  \item{pkgPath}{\code{pkgPath} a character string for the path where
    the data package created will be stored}
  \item{version}{\code{version} a character string for the version
    number of the data package to be created}
  \item{author}{\code{author} a list with an author elementfor the name
    of the author of the data package and a maintainer element for the
    name and email address of the maintainer of the dat package to be
    created} 
  \item{fromWeb}{\code{fromWeb} a boolean indicating whether the data
    will be read from the internet or locally}
  \item{url}{\code{url} an URL of file name to read the data from}
  \item{key}{\code{key} a character string for the name of Swiss-Prot
    annotation element, e. g. "Swiss-Prot accession number"}
  \item{envName}{\code{envName} a character string for the name of an
    environment object} 
}
\details{
  If \code{fromWeb} is FALSE, url will be the file name of a local file.
}
\value{
  This function returns NULL
}
\references{\url{ftp://ftp.ebi.ac.uk/pub/databases/swissprot/release/sprot41.dat}} 
\author{Jianhua Zhang}

\seealso{\code{\link{ABPkgBuilder}}}
\keyword{manip}


