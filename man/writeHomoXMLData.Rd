\name{writeHomoXMLData}
\alias{writeHomoXMLData}
\alias{writeHGID2Caption}
\alias{writeHGID2LL}
\alias{writeHomoData}
\alias{homoXMLParser}
\title{Functions to parse HomoloGene XML data file and build the
  homology annotation data package}
\description{
  HomoloGene maintains a homology XML data file that differes both in
  the format and contents from the old text file version. The functions
  described here parse the file build the homology annotation data
  pacakge based on the source data.  
}
\usage{
writeHomoXMLData(pkgName = "homology", pkgPath, version, author, url =
"ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build39.2/homologene.xml.gz")
writeHGID2Caption(pkgName, pkgPath, hgid2Cap)
writeHGID2LL(pkgName, pkgPath, hgid2LL)
writeHomoData(pkgName, pkgPath, homoFile)
homoXMLParser(fileName)
}

\arguments{
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be created}
  \item{pkgPath}{\code{pkgPath} a character string for the path to the
    directory where the new data package will be stored}
  \item{version}{\code{version} a character string for the version
    number of the data package to be created}
  \item{author}{\code{author} a list of character strings with an author
    (name of the author of the data package) and a maintainer (e-mail of
    the author of the package)}
  \item{url}{\code{url} a character string for the url to the ftp site
    where the HomoloGene XML file is available. The url change with
    different builds. Check the HomoloGene web site for the latest one}
  \item{hgid2Cap}{\code{hgid2Cap} a matrix containing mappings between
    HGIDs and their textual descriptions}
  \item{hgid2LL}{\code{hgid2LL} a matrix with mappings between HGIDs and
    LocusLink ids}
  \item{homoFile}{\code{homoFile} a character string for the name of a
    file containing data for homologous genes}
  \item{fileName}{\code{fileName} a character string for the name of the
    XML file downloaded/unzipped from HomoloGene's ftp site}
}
\details{
  \code{\link{writeHomoXMLData}} calls other functions listed in this
  help page to complete it's tasks. All the other functions are help
  functions that may not of greate interest to users.
}
\value{
  The function returns invisible(NA)
}
\references{\url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=homologene}}
\author{Jianhua Zhang}

\keyword{manip}
