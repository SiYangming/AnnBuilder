\name{writeManPage}
\alias{writeManPage}
\alias{writeMan4Fun}
\alias{formatName}
\alias{writeREADME}
\alias{writeDescription}
\alias{getDsrc}
\alias{getItem}
\alias{getSrcNBuilt}
\alias{getUrlNBuilt}
\alias{getDSrc}
\alias{writeAccessory}
\alias{writeZZZ}
\alias{getAllRdaName}
\alias{writeFun}
\alias{escapeLatexChr}
\alias{writeMan4QC}
\alias{getExample}
\alias{getSrcBuiltNRef}
\alias{getBuild4Yeast}
\title{Functions that write supporting files needed by a data package}
\description{
  The functions are mainly used to write man pages and supporting
  functions that are needed for a data package  
}
\usage{
writeManPage(pkgName, pkgPath, manName, organism = "human", src = "ll",
isEnv = TRUE) 
writeMan4Fun(pkgName, pkgPath, organism = "human", QCList, dSrc = "all" )
formatName(toFormat)
writeREADME(pkgPath, pkgName, urls)
writeDescription(pkgName, pkgPath, version, author, dataSrc, license)
getDSrc(organism)
getSrcNBuilt(dSrc, organism)
getUrlNBuilt(src, organism)
getDSrc{organism}
writeAccessory(pkgName, pkgPath, organism, version, author = list(author =
"who", maintainer = "who@email.net"), dataSrc, license)
writeFun(pkgPath, pkgName, organism = "human")
writeZZZ(pkgPath, pkgName)
getAllRdaName(pkgName, pkgPath)
escapeLatexChr(item)
writeMan4QC(pkgName, pkgPath)
getExample(pkgName, manName, isEnv = TRUE)
getSrcBuiltNRef(src, organism)
getBuild4Yeast(src, manName)
}

\arguments{
  \item{pkgName}{A character string for the name of a data package or R
    library}
  \item{pkgPath}{A character string for the path where pkgname resides}
  \item{organism}{A character string for the name of the organism of
    interests} 
  \item{toFormat}{A character string form whom any underscore will be
    removed}
  \item{urls}{A vector of character of string for the urls of the data
    source used to create the rda files}   
  \item{dSrc}{A vector of character strings containing the short names
    of public data sources (e. g. LL for LocusLink)}
  \item{src}{A character string for the short name of a public data
    source}
  \item{version}{A character string for the version number}
  \item{author}{A named vector of strings with two elements named name
    and address, respectively. Name is a character string for the name
    of the person who maintains the data package and address is the email
    address of the person}
  \item{item}{A character string to be escaped by if it is a latex
    character}
  \item{QCList}{A list with statistical data derived from
    \code{\link{getDPStats}}}
  \item{manName}{\code{manName} a character string for the name of the
    man page to be created}
  \item{isEnv}{\code{isEnv} a boolean to indicate whether the object a
    man page concerns is an R environment or not}
  \item{dataSrc}{\code{dataSrc} a vector of character strings for the
    data sources used to create a package}
  \item{license}{\code{license} a character string for the license the
    package is under}
}
\details{
  If pkgname = "XX" and elenames = "yy", the Rd file will be "XXyy.Rd"
  appended to the path if short is FALSE. Otherwise, the Rd file will be
  "yy.Rd" appended to the path.

  \code{\link{writeManPage}} writes a man page for a given object that
  is stored in the data directory.

  \code{\link{getExample}} creates a set of example code that is going
  to be used in a man page depending on whether the man page is for an
  environment object or not.
  
  \code{\link{getSrcBuiltNRef}} creates the text that is going to be
  used for built and reference information in a man page.

  \code{\link{getBuild4Yeast}} creates the text that is going to be used
  for built and reference information for the man page for yeast.
}
\value{
  All functions return a character string.
}
\references{An Introduction to R - Writing R Extensions}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide BioInformatics functionalities through R}

\seealso{\code{\link{ABPkgBuilder}}} 

\examples{
\dontrun{
makeSrcInfo()
dir.create(file.path(".", "pkg"))
dir.create(file.path(".", "pkg", "data"))
dir.create(file.path(".", "pkg", "man"))
writeManPage("pkg", getwd(), "CHR")
list.files(file.path(getwd(), "pkg", "data"))
unlink("pkg", TRUE)
}
}
\keyword{manip}
