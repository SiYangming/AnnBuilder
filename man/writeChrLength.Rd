\name{writeChrLength}
\alias{writeChrLength}
\alias{findChrLength}
\alias{writeOrganism}
\title{Functions that creates binary files for chromosome length and organism}
\description{
  These functions figure out the chromosome length and write the length
  and organism binary files to the data directory of the pacakge
}
\usage{
writeChrLength(pkgName, pkgPath, chrLengths)
findChrLength(organism, srcUrl = getSrcUrl("GP", organism))
writeOrganism(pkgName, pkgPath, organism)
}

\arguments{
  \item{pkgName}{\code{pkgName} a character string for the name of a
    data package or R library}
  \item{pkgPath}{\code{pkgPath} a character string for the path where
    pkgname resides} 
  \item{organism}{\code{organism} a character string for the name of the
    organism of interests}
  \item{srcUrl}{\code{srcUrl} a character string for the url of the data
    source used to create the binary file for chromosome length}   
  \item{chrLengths}{\code{chrLengths} a named vector of integers with
    the names being the chromosome numbers and the values of the vector
    being the total lengths of chromosomes}
}
\details{
  \code{\link{findChrLength}} extracts data from the source and figures
  out the total length for each chromosome. The total length for a
  chromosome is determined as the maximum chromosome location plus 1000.

  \code{\link{writeChrLength}} writes the chromosome length data to the
  data directory as a binary file.

  \code{\link{writeOrganism}} writes the name of the organism to the
  data directory as a binary file.
}
\value{
  \code{\link{findChrLength}} returns a named vector of integers. 
}
\references{}
\author{Jianhua Zhang}

\seealso{\code{\link{ABPkgBuilder}}} 

\examples{
\dontrun{
    path <- file.path(.path.package("AnnBuilder", "temp"))
    dir.create(file.path(path, "test"))
    dir.create(file.path(path, "test", "data"))
    chrLength <- findChrLength("human")
    writeChrLength("test", path, chrLength)
    writeOrganism("test", path, "human")
    list.files(file.path(path, "test", "data"))
    unlink(file.path(path, "test"), TRUE)
}
}
\keyword{manip}
