\name{homoPkgBuilder}
\alias{homoPkgBuilder}
\alias{procHomoData}
\alias{getLL2IntID}
\alias{getIntIDMapping}
\alias{mapIntID}
\alias{mapOrgs}
\alias{writeRdaNMan}
\alias{mapPS}
\alias{getHomoPS}
\alias{getHomoDList}
\alias{getHomoData}
\alias{getHomoData}
\alias{getOrgNameNCode}
\alias{saveOrgNameNCode}
\alias{HomoData2List}
\alias{orgNameNCode}
\title{Functions to build a homology data package using data from NCBI}
\description{
  This function builds a data package that maps internal HomoloGene ids
  of an organism to LocusLink ids, UniGene ids, percent identity of the
  alignment, type of similarities, and url to the source of a curated
  orthology of organisms of all pairwise best matches based on data from
  \url{ftp://ftp.ncbi.nih.gov/pub/HomoloGene/hmlg.ftp}
}
\usage{
homoPkgBuilder(suffix = "homology", pkgPath, version, author, url =
getSrcUrl("HG"))
procHomoData(url = getSrcUrl("HG"))
getLL2IntID(homoData, organism = "")
mapOrgs(toMap, what = c("code","name"))
mapPS(homoMappings, pkgName, pkgPath, tempList)
getHomoDList(data, what = "old")
getHomoData(entries, what = "old", objOK = FALSE)
getOrgNameNCode()
saveOrgNameNCode(pkgName, pkgPath, tepList)
HomoData2List(data, what = "old")
}
\arguments{
  \item{suffix}{\code{suffix} a character string for the suffix to be
    attached to the end of a three-letter short form for an organism to
    form the name of a package to be created for homologous genes of the
    organism}
  \item{pkgName}{\code{pkgName} a character string for the name of data
    package to be built}
  \item{pkgPath}{\code{pkgPath} a character string for the name of the
    directory where the created package will be stored}
  \item{version}{\code{version} a character string for the verion number
    of the package to be built}
  \item{author}{\code{author} a list with an author element for the name
    of the author and a maintainer element for the name and e-mail
    address of the maintainer of the package} 
  \item{url}{\code{url} the url to the ftp site from which the source data
    file can be obtained. The default value is
    url{ftp://ftp.ncbi.nih.gov/pub/HomoloGene/hmlg.ftp}}
  \item{homoData}{\code{homoDatgetHomoDLista} a data frame that contains
    the homology data from the source}
  \item{homoMappings}{\code{homoMappings} same as homoData but only
    contains data for an organism of concern}
  \item{organism}{\code{organism} a character string for the name of the
    organism of interest}
  \item{toMap}{\code{vect} a vector of character strings}
  \item{entries}{\code{entries} a vector of character strings}
  \item{data}{\code{data} a data matrix}
  \item{what}{\code{what} a character string that can either be "code"
    or "name" for function mapOrgs or "old" or "xml" for functions
    getHomoDList, getHomoData, and HomoData2List}
  \item{tepList}{\code{tepList} a list containing key and value pairs
    that are going to be used to replace the corresponding matching
    items in a template file for man pages}
  \item{tempList}{\code{tempList} same as \code{tepList}}
  \item{objOK}{\code{objOK} a boolean indicating whether the homoDATA
    environment will be a list of homoDATA (TRUE) objects or lists (FALSE)}
}
\details{
  procHomoData process the source data and put the data into a data
  frame that will be used later.

  getLL2IntID maps LocusLink ids to HomoloGene internal ids

  getIntIDMapping maps HomoloGene ids to ids include LocusLink ids,
  GneBank accession numbers, percent similarity values, type of
  similarities, and the url to the curated orthology.

  mapIntID captures the reverse mapping between reciprocal homologous
  genes.

  mapOrgs converts organism codes to scientific names.

  writeRdaNMan creates an rda file and the corresponding man page for a
  data environment.

  mapPS maps HomologGene Internal ids to homoPS objects generated using
  data from the source.

  getHomoPS creates a homoPS object using data passed as a vector.
}
\value{
  procHomoData, mapIntID, and getLL2IntID returns a matrix.

  getIntIDMapping returns an R environment with mappings between
  HomoloGene internal ids and mapped data.

  getHomoPS returns a homoPS object with slots filled with data passed.

  mapOrgs returns a vector of character strings.
}
\references{\url{ftp://ftp.ncbi.nih.gov/pub/HomoloGene/README}}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities using R}

\seealso{\code{\link{ABPkgBuilder}}}
\keyword{manip}
