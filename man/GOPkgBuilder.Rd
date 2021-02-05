\name{GOPkgBuilder}
\alias{GOPkgBuilder}
\alias{writeDocs}
\alias{copyTemplates}
\alias{getRepList}
\title{Functions to build a data package using GO data}
\description{
  
  WARNING: DO NOT use this function. Use Bioconductor packages such as
  biomaRt to obtain GO mappings more recent than those available in the
  current Bioconductor release.

  \url{http://www.bioconductor.org/}
 
  These functions creates data, documentation, and other supporting
  files that consist an annotation data package using data from GO. 
}
\usage{
GOPkgBuilder(pkgName, pkgPath, filename, version, author, lazyLoad=TRUE)
writeDocs(baseName, pkgName, pkgPath, version, author, repList, pattern,
isFile = TRUE) 
copyTemplates(repList, pattern, pkgName, pkgPath, replaceBy = NULL)
getRepList(what, srcObjs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be built}
  \item{pkgPath}{Describe \code{pkgPath} a character string for the path
    to which the data package to be built will be stored}
  \item{filename}{Name of the GO file to parse.  This file should be
  from the GO website in OBO XML format.}
  \item{version}{\code{version} a character string for the version
    number of the data package}
  \item{author}{\code{author} a named vector of character string with a
    name element for the name of the author and address element for the
    email address of the author}
  \item{repList}{\code{repList} a list with LLSOURCE, GOSOURCE,
    LLBUILT, GOBUILT, and DATE elements containing source url or built
    date information that will be used to replace corresponding texts in
    man page templates stored in the templates subdirectory}
  \item{pattern}{\code{pattern} a character string that will
    be used as a pattern to copy man page templates files in the "templates"
    subdirectory to the "man" subdirectory of a newly created data
    package using the function copySubstitute of Biobase}
  \item{replaceBy}{\code{replaceBy} a character string specifying the
    text used to replace the pattern contained by the name of a template
    man page files when writing to a newly created data package}
  \item{what}{\code{what} a character string for the name of the data
    package to be created for which a replacement list will be
    generated}
  \item{baseName}{\code{baseName} a character string for the name of the
    base file based on which a data package was built. "" if there is
    none}
  \item{srcObjs}{\code{srcObjs} a list containing source data objects
    that are sub classes of pubRepo}
  \item{isFile}{\code{isFile} a boolean indicating whether
    \code{baseName} is a file or an R object}
  \item{lazyLoad}{\code{lazyLoad} a boolean indicating whether a lazy
    load database will be created}
}
\details{
  This package relies on the xml data file from
  \url{http://www.godatabase.org/dev/database/archive/2003-04-01/go_200304-termdb.xml.gz} to obtain the data. The url changes when data are updated. The system has built in code to figure out where the latest data are and use that data to build the data package. 
}
\value{
  This function does not return any value
}
\references{\url{http://www.godatabase.org}}
\author{Jianhua Zhang}

\seealso{\code{\link{ABPkgBuilder}}, \code{\link{KEGGPkgBuilder}}}
\examples{
if(interactive()){
GOPkgBuilder(pkgName = "GO", pkgPath = tempdir(), version = "1.2.1",
goUrl = "http://www.bioconductor.org/datafiles/wwwsources/Tgo.xml",
author = c(author = "who", maintainer = "who@email.com"))
list.files(file.path(tempdir(), "GO"))
unlink(file.path(tempdir(), "GO"), TRUE)
}
}
\keyword{manip}

