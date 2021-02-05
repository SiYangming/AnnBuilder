\name{getDPStats}
\alias{getDPStats}
\alias{getDate}
\alias{getProbeNum}
\alias{matchProbes}
\alias{getPBased}
\alias{formatABQCList}
\alias{countMapping}

\title{Functions to read in the statistics about a data package}
\description{
  These functions generate a list showing the name, data of creation,
  number of genes for each rda file, and the actual number of genes that
  get mapped for each rda file.
}
\usage{
getDPStats(baseF, pkgName, pkgPath, saveList = TRUE, isFile = TRUE)
getDate(pkgName, pkgPath, fromDesc)
getProbeNum(pkgName, pkgPath, noNA = FALSE)
matchProbes(baseF, pkgName, pkgPath, toMatch, isFile = TRUE)
getPBased()
formatABQCList(x)
countMapping(rdaName, noNA = FALSE)
}

\arguments{
  \item{baseF}{\code{baseF} a character string for the name of a file
    that is going to be used as the base file to calculate the total
    number of probes and matched probes by a data package. Set to "" if
    there is no base file}
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package of concern} 
  \item{pkgPath}{\code{pkgPath} a character string for name of the path
    to which the data package is stored.}
  \item{noNA}{\code{noNA} a boolean to indicate whether counts will
    exclude entries with NA as the value.}
  \item{saveList}{\code{sageList} a boolean indicating whether the
    results will be returned as a list only (FALSE) or saved to a file as
    well (TRUE)} 
  \item{toMatch}{\code{toMatch} a vector of character strings for the
    names of the rda files whose keys will be matched against the probe
    ids of a base file (baseF)}
  \item{x}{\code{x} a list object produced by function
    \code{\link{getDPStats}}} 
  \item{fromDesc}{\code{fromDesc} a boolean that will get a date from a
    DESCRIPTION file if set TRUE or the current date if FALSE}
  \item{isFile}{\code{isFile} a boolean that will be TRUE if
    \code{baseF} is the name of a file}
  \item{rdaName}{\code{rdaName} a character string for the name of an
    rda file whose man page will be generated}
}
\details{
  Date of creation is the date when the package was created using
  AnnBuilder and in most cases is not the date when the source file
  AnnBuilder used to create the rda files was created. The date when the
  source data were built are listed in the man page for the package
  (?package name).

  The number of genes and number of genes mapped normally differ because
  not all genes in a given set can be mapped to annotation data. For
  probe based rda files (e. g. maps Affymetrix ids to annotation data), the
  number of mapped genes out of the total is given. For non-probe based
  rda files, only the total number of mapped items is given.

  The total number of probes of each rda file will be checked against
  the total of the base file and the names of the rda files whose total
  is off will be listed.
}
\value{
  \item{list}{A list with name and value pairs}
}
\author{Jianhua Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{ABPkgBuilder}}}

\examples{
# Run this code after changing the settings correctly
# Change the varaibles before you run the code
pkgName <- "hgu95a"
pkgPath <- "where/your/data/package/is"
# Call getABStats
# getDPStats(pkgName, pkgPath)
}
\keyword{misc}




