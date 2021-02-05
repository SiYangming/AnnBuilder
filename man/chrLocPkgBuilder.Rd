\name{chrLocPkgBuilder}
\alias{chrLocPkgBuilder}
\alias{getChrNum}
\alias{saveCytoband}
\alias{getChroms4Org}
\title{A function to build a data package containing mappings between
  LocusLink ids and the chromosomal locations of genes represented by
  the LocusLink ids}
\description{
  This function uses data provided by UCSC to build a data package that
  contains mappings between LocusLink ids and chromosome numbers and the
  chromosomal location of genes represented by LocusLink ids on each
  chromosome
}
\usage{
chrLocPkgBuilder(pkgName = "humanCHRLOC", pkgPath, version, author,
organism = "Homo sapiens")
getChrNum(chr)
saveCytoband(pkgName, pkgPath, organism, url, ext = "cytoBand.txt.gz")
getChroms4Org(organism)
}
\arguments{
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be created}
  \item{pkgPath}{\code{pkgPath} a character string for the directory
    where the created data package will be stored}
  \item{version}{\code{version} a character string for the version
    number of the data package to be created}
  \item{author}{\code{author} a list with an author element for the name
    of the creater of the data package and a maintainer element for the
    email address of the creater} 
  \item{organism}{\code{organism} a character string for the organism
    of concern}
  \item{url}{\code{url} a character string of the url of UCSC ftp site
    where to file refLink.txt.gz and refGene.txt.gz are stored. The
    files will be used to produce the data package}
  \item{chr}{\code{chr} a character string for the chromosome number
    extracted from the source data}
  \item{ext}{\code{ext} a character string for the file name to be
    appanded to the argument url}
}
\details{
  The data package created maps LocusLink ids to chromosomal
  locations. Mappings of other public data repository ids including Gene
  Ontology, RefSeq, and UniGene to LocusLink ids can be made available
  using \code{\link{map2LL}} 
}
\value{
  invisible
}
\references{}
\author{Jianhua Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institue to provide bioinfomatics functionalities through R}
\seealso{\code{\link{map2LL}}}
\examples{
# Please note that the example will take a while to finish
if(interactive()){
chrLocPkgBuilder(pkgName = "humanCHRLOC", pkgPath = tempdir(),
version = "1.0.1", author = list(author = "who", maintainer =
"who@email.com"), organism = "human")
}
}
\keyword{manip}
