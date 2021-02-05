\name{pfamBuilder}
\alias{pfamBuilder}
\title{Building Functions for the Data Package of Pfam Database}
\description{ 
  These functions builds a data package for Pfam database
}
\usage{
pfamBuilder(pkgName="PFAM", pkgPath, version, author,
            fromWeb=TRUE, lazyLoad=TRUE, useTmp=FALSE, sqlFile=NULL)
}
\arguments{
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be built.  The default is "PFAM"}
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
  \item{useTmp}{\code{useTmp} a boolean.  If TREU, the sqlite file will
    be saved to tempdir.  Otherwise, it will be saved to the data
    subdirectory of the package.  The default is FALSE.}
  \item{sqlFile}{\code{sqlFile} a character string to indicate the full
    path of the local SQL file.  The default is NULL so that the
    SQL file will be downloaded from the Pfam website}

}
\references{\url{http://www.sanger.ac.uk/Software/Pfam/}}
\author{Ting-Yuan Liu}
\examples{
# Not provided. 
}
\keyword{manip}
