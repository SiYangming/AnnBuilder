\name{map2LL}
\alias{map2LL}
\alias{getExten}
\alias{saveColSepData}
\alias{getOrgName}
\alias{getReverseMapping}
\alias{saveData2Env}
\alias{reverseMap4GO}
\alias{getLL2ACC}

\title{A function that maps LocusLink ids to other public repository ids
and vice versa} 
\description{
  This function uses data files provided by NCBI to create a data
  package that contains mappings between LocusLink ids and GO, RefSeq,
  and UniGene ids and vice versa 
}
\usage{
map2LL(pkgName, pkgPath, organism, version, author, eg = EG(parser =
file.path(.path.package("AnnBuilder"), "scripts", "egLLMappingUGParser")), lazyload = TRUE)
getExten(what)
getOrgName(organism, what = c("common", "scientific"))
getReverseMapping(data, sep = ";")
saveData2Env(data, fun = splitEntry, pkgName, pkgPath, envName)
reverseMap4GO(data, sep = ";", type = c("ll2GO", "GO2LL")  )
getLL2ACC(url = paste("ftp://ftp.ncbi.nih.gov/refseq/LocusLink/",
getExten("acc"), sep = ""), organism = "human") 
}
\arguments{
  \item{organism}{\code{organism} a character string for the name of the
    organism of interest}
  \item{pkgPath}{\code{pkgPath} a character string for the name of the
    directory where the created data package will be stored}
  \item{version}{\code{version} a character string for the version
    number of the data package to be created}
  \item{author}{\code{author} a list with an author element for the
    name of the creater of the data package and a maintainer element for
    the email address of the creater} 
  \item{url}{\code{url} a character string for the url of NCBI's ftp
    site where source data are stored. Current value is
    \url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/}}
  \item{what}{\code{what} a character string for the type of mapping
    source data (i. e. "go", "ug" ...) or description of organism
    name("scientific" or "short")}
  \item{data}{\code{data} a matrix to be processed}
  \item{sep}{\code{sep} a character string the separator used to
    separate data elements for a given entry}
  \item{envName}{\code{envName} a character string for the name of the
    environment object to be stored in the data package to be created}
  \item{fun}{\code{fun} the name of an R function to be called to
    process a data set before storing the data to an environment object}
  \item{pkgName}{\code{pkgName} a character string for the name of data
    package to be created}
  \item{type}{\code{what} a character string that should either be
    "ll2GO" or "GO2LL" to indicate a reverse mapping from LocusLink id
    to GO or vice versa}
  \item{lazyload}{\code{lazyload} a boolean indicating whether a lazy
    load database will be created}
  \item{eg}{\code{eg} an EG object}
}
\details{
  Three files namely loc2go, loc2ref, and loc2UG will be used to create
  the mappings. The files were in
  \url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/} at the time of the
  writing. \code{\link{getExten}} maintains names for the three
  files. Should any of the names been changed by the server,
  \code{\link{getExten}} has to be modified.

  \code{\link{getExten}} and \code{\link{saveColSepData}} are supporting
  functions to \code{\link{map2LL}} 
}
\value{
  invisible
}
\references{\url{http://www.ncbi.nlm.nih.gov/LocusLink/}}
\author{Jianhua Zhang}

\examples{
\dontrun{
# Please note that the example will take a while to finish
map2LL(pkgPath = tempdir(), version = "1.0.0", organism = "human",
       author = list(author = "myName", maintainer = "myeail@email.com"))
}
}
\keyword{manip}



