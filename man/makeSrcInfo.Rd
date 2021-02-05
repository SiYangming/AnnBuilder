\name{makeSrcInfo}
\alias{makeSrcInfo}
\alias{getAllSrc}

\title{Functions to make source information available for later use}
\description{
  These functions read from a text file (AnnInfo) that have been stored
  in the data directory and create an environment object called AnnInfo
  that will be available for later access 
}
\usage{
makeSrcInfo(srcFile = "")
getAllSrc()
}

\arguments{
  \item{srcFile}{\code{srcFile} a character string for the name of the
    source file that contains source data information}
}
\details{
  The environment object created (AnnInfo) is a list with four elements:
  \item{short}{a character string for the description that will be used
    to describe an annotation element in an XML file to be generated}
  \item{long}{a character string that will be used to describe an
    annotation element in the help file for a given data environment
    that will be contained in a data package to be created}
  \item{src}{a character string for the short hand name of the source
    (e.g. ll for LocusLink)}
  \item{pbased}{a boolean that is TRUE if the annotation element is for
    a probe or FALSE otherwise}
}
\value{
  \code{\link{getAllSrc}} return a vector of character string for short
  hand names of data sources
}
\author{Jianhua Zhang}

\seealso{\code{\link{ABPkgBuilder}}, \code{\link{GOPkgBuilder}},
  \code{\link{KEGGPkgBuilder}}} 
\examples{
\dontrun{
  makeSrcInfo()
  ls(AnnInfo)
}
}
\keyword{file}

