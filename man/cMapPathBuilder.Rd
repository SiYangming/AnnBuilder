\name{cMapPathBuilder}
\alias{cMapPathBuilder}
\alias{cMAPParser}
\title{Functions that build a data pacakge using data provided by cMAP}
\description{
  NCICB Pathway Interaction Database provides two data files molecule
  interaction data for BioCarta and KEGG pathways. The functions
  described here use the two files and build a data package containing
  the data
}
\usage{
cMapPathBuilder(cartaName, keggName, pkgName = "cMAP", pkgPath, version
= "1.1.0", author = list(author = "anonymous", maintainer =
"anonymous@email.com"), lazyLoad = TRUE)
cMAPParser(sourceFile)  
}

\arguments{
  \item{cartaName}{\code{cartaName} a character string for the name of
    the XML file containing data for BioCarta pathways}
  \item{keggName}{\code{keggName} a character string for the name of
    the XML file containing data for BioCarta pathways}
  \item{pkgName}{\code{pkgName} a character string for the name of the
    package to be built}
  \item{pkgPath}{\code{pkgPath} a character string for the path to the
    directory where the new package to built will reside}
  \item{version}{\code{version} a character string for the version
    number of package to be built}
  \item{author}{\code{author} a list with an author and maintainer
    element for the name and email address of the author of the pacakge}
  \item{sourceFile}{\code{sourceFile} a character string for the name of
    the source data for BioCarta or KEGG}
  \item{lazyLoad}{\code{lazyLoad} a boolean indicating whether a lazy
    load database will be created}
}
\details{
  cMAP currently does not support ftp downloading of the source data
  file. The files to be used have to be downloaded through a web browser
  and the name (with full path) of the downloaded file will be used.
}
\value{
  cMapPathBuilder returns a list with three elements:
  \item{molecule}{a list of vectors/lists containing molecule data}
  \item{interaction}{a list of vectors/lists containing molecule
    interaction data}
  \item{pathway}{a list of vectors containing pathway component data}
}
\references{\url{http://cmap.nci.nih.gov/PW/}}
\author{Jianhua Zhang}
\note{}

\examples{
## No example is provided
}
\keyword{manip}
