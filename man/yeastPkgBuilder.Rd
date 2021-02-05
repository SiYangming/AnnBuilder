\name{yeastPkgBuilder}
\alias{yeastPkgBuilder}
\alias{findYGChrLength}
\title{Functions to do a data package for yeast genome}
\description{ 
  These functions builds a data package for yeast genome using data from
  Yeast Genome web site of Stanford University, KEGG, and Gene Ontology.
}
\usage{
yeastPkgBuilder(pkgName, pkgPath, base="", version, author,
                            fromWeb=TRUE, lazyLoad=TRUE)
}
\arguments{
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be built}
  \item{base}{\code{base} a matrix with two columns with the first one
    being probe ids and the second one being their mappings to ORF (Open
    Reading Frame) ids. Columns have the name "probe" and "orf"}
  \item{pkgPath}{\code{pkgPath} a character string for the directory
    where the data package to be built will be stored}
  \item{version}{\code{version} a character string for the version
    number of the data package to be built}
  \item{author}{\code{author} a named vector of two character strings
    with a name element for the name and an address element of email
    address of the maintainer of the data package}
  \item{fromWeb}{\code{fromWeb} a boolean to indicate whether the data
    from GO should be downloaded from the web or read locally. The url
    for GO should be the file name of a local file if fromWeb is
    FALSE. For windows users, the data file from GO should be
    downloaded/unzipped manually and set the url for GO to be the name
    of the local file}
  \item{lazyLoad}{\code{lazyLoad} a boolean indicating whether a lazy
    load database will be created}

}
\details{
  Annotation elements are limited to those provided by Yeast Genome
  (gene name, chromosome number, chromosomal location, GO id, and
  evidence code), KEGG (path and enzyme data) and GO (GO mappings)
}
\value{
  \code{findYGChrLength} returns a named vector of integers with
  chromosome numbers names and length of chromosomes as values. 
}
\references{\url{http://www.yeastgenome.org}}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}
\examples{
# Not provided. 
}
\keyword{manip}
