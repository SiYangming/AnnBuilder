\name{KEGGPkgbuilder}
\alias{KEGGPkgBuilder}
\alias{getEIdNName}
\alias{getKEGGFile}
\alias{getKEGGGeneMap}

\title{A function to make the data package for KEGG}
\description{
  This function generates a data package with rda files mapping KEGG
  pathway or enzyme names to ids and vice versa. The source files for
  making the mapping are from the Internet.
}
\usage{
KEGGPkgBuilder(pkgPath, pkgName = "KEGG", version = "1.0.1", 
author = list(author = "who",maintainer = "who@email.com"))
getEIdNName(enzymeURL)
getKEGGFile(whichOne, organism = "hsa")
getKEGGGeneMap(organism = "Homo sapiens"))
}

\arguments{
  \item{pkgPath}{A character string for the name of path to which the
    data package will be stored.}
  \item{pkgName}{A character string for the name of the data package.}
  \item{version}{A character string for the version number of the system
    by which the data package is generated.}
  \item{author}{A list of character strings with one element being name
    for the name of the author and another being address being the email
    address of the author}
  \item{organism}{\code{organism} a character string for the name of the
    organism of interest}
  \item{whichOne}{A charcter string for the name of file type. Valid
    values include "path" or "enzyme"}
  \item{enzymeURL}{A character string for the URL from which the source
    file for enzyme data will be downloaded.}
}
\details{
  The data package produced will have the normal structure of an R
  package (i. g. with R, man, data, and src directories) under a
  directory defined by pkgName under pkgPath.
}
\value{
  This function does not return any value.
}
\references{An Introduction to R - Writting R Extensions.}
\author{Jianhua Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R.}

\seealso{\code{\link{package.skeleton}}}

\keyword{manip}

