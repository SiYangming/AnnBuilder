\name{ABPkgBuilder}
\alias{ABPkgBuilder}
\alias{getBaseParsers}
\alias{createEmptyDPkg}
\alias{getDirContent}
\alias{getMultiColNames}
\alias{getUniColNames}
\alias{getTypeColNames}
\alias{splitEntry}
\alias{twoStepSplit}
\alias{saveMat}
\alias{saveList}
\alias{getChrLenghts}
\alias{nameGOByCat}
\alias{getChrLengths}
\alias{getHumanChrLengths}
\alias{getMouseChrLengths}
\alias{getRatChrLengths}
\alias{getYeastChrLengths}
\alias{getList4GO}
\alias{vect2List}
\alias{resumeSrcUrl}
\alias{writeDatalist}
\alias{getEGAccName}

\title{Functions that support a single API for building data packages}
\description{
  These functions support a single API represented by ABPkgBuilder to
  allow users to build annotation data packages by providing a limited
  number of parameters. Other parameters will be figured out by the
  supporting functions. 
}
\usage{
ABPkgBuilder(baseName, srcUrls, baseMapType = c("gb", "ug", "ll", "image",
"refseq", "gbNRef"), otherSrc = NULL, pkgName, pkgPath, organism,
version, author, fromWeb = TRUE, lazyLoad = TRUE)
getBaseParsers(baseMapType = c("gb", "ug", "image", "ll", "refseq",
"gbNRef", "ll2gb", "gb2ll", "eggo", "eginfo", "egrefseq", "egpubmed", "egunigene", "egmim")))
createEmptyDPkg(pkgName, pkgPath, folders, force = TRUE)
getDirContent(dirName, exclude = NULL)
getMultiColNames()
getUniColNames()
getTypeColNames()
splitEntry(dataRow, sep = ";", asNumeric = FALSE)
twoStepSplit(dataRow, entrySep = ";", eleSep = "@", asNumeric = FALSE)
saveMat(data, pkgName, pkgPath, envName, keyCol = 1,
                         valCol = 2, fun = function(x) x)
saveList(dList, pkgName, pkgPath, envName)
nameGOByCat(GOWithEvi, goCat)
getChrLengths(organism)
getHumanChrLengths()
getMouseChrLengths()
getRatChrLengths()
getYeastChrLengths()
getList4GO(goNCat, goNEvi)
vect2List(vector, vectNames)
resumeSrcUrl(srcObjs, organism)
writeDatalist(pkgName, pkgPath)
getEGAccName()
}

\arguments{
  \item{baseName}{\code{baseName} a character string for the name of a
    file to be used as a base file to base source data. The file is
    assumed to have two columns (separated by tabs "\t") with the first
    one being the names of genes (probes) to be annotated and the second one
    being the maps to GenBank accession numbers, UniGene ids, image
    clone ids or LocusLink ids} 
  \item{srcUrls}{\code{srcUrls} a vector of named character strings for
    the urls where source data files will be retrieved. Valid sources are
    LocusLink, UniGene, Golden Path, Gene Ontology, and KEGG. The names
    for the character strings should be LL, UG, GP, GO, and KEGG,
    respectively. LL and UG are required} 
  \item{baseMapType}{\code{baseMapType} a character string that is
    either "gb","ug", "image", "ll", "image", "refseq", "gbNRef" to
    indicate whether the probe ids in baseName are mapped to GenBack
    accession numbers, UniGene ids, image clone ids, LocusLink ids,
    RefSeq ids, or a mixture of GenBank accession numbers and RefSeq ids}
  \item{otherSrc}{\code{otherSrc} a vector of named character strings
    for the names of files that contain mappings between probe ids of
    baseName and LobusLink ids that will be used to obtain the unified
    mappings between probe ids of baseName and LocusLink ids based on
    all the sources. The strings should not contain any number and the
    files have the same structure as baseName}
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be built (e. g. hgu95a, rgu34a)}
  \item{pkgPath}{\code{pkgPath} a character string for the full path of
    an existing directory where the built package will be stored}
  \item{organism}{\code{organism} a character string for the name of the
    organism of concern (now can only be "human", "mouse", or "rat")}
  \item{version}{\code{version} a character string for the version number}
  \item{author}{\code{author} a list of character strings with an author
    element for the name of the author and maintainer element for the
    email address of the author}
  \item{force}{\code{force} a boolean that is set to TRUE if the package
    to be created will replace an existing package with the same name}
  \item{dirName}{\code{dirName} a character string for the name of a
    directory whose contents are of interests}
  \item{exclude}{\code{exclude} a character string for a pattern matching
    parameter that will be used to exclude contents of a directory that
    match the pattern}
  \item{dataRow}{\code{dataRow} a character string containing data
    elements with elements separated by \code{sep} or \code{entrySep}
    and a descriptive string attached to each element following
    \code{eleSep}}
  \item{sep}{\code{sep} a character string for a separator}
  \item{entrySep}{\code{entrySep} a character string for a separator}
  \item{eleSep}{\code{eleSep} a character string for a separator}
  \item{asNumeric}{\code{asNumeric} a boolean that is TRUE when the
    splited values will be returned as numeric values}
  \item{fromWeb}{\code{fromWeb} a boolean to indicate whether the source
    data will be downloaded from the web or read from a local file}
  \item{folders}{\code{folders} a vector of character strings for the
    names of folders to be created within a package that is going to be
    created}
  \item{data}{\code{data} a data matrix to be written as an environment
    object}
  \item{dList}{\code{dList} a list to be written an an environment
    object}
  \item{envName}{\code{envName} a character string for the name of an
    environment object to be written as keys in an environment}
  \item{keyCol}{\code{keyCol} a numeric number indicating the column of
    a matrix that contains keys}
  \item{valCol}{\code{valCol} a numeric number indicating the column of
    a matris that contains data that will be written as values in an
    environmnet}
  \item{fun}{\code{fun}an R function that will be passed as an argumnet}
  \item{GOWithEvi}{\code{goWithEvi} a vector of character string in the
    format of "GO:xxxx@TS;GO:xxxxx@P;..." where letters following "@"
    are evidence code}
  \item{goCat}{\code{goCat} a matrix with the first column being GO ids
    and the second column being GO categories}
  \item{goNCat}{\code{goNCat} a named vector with GO category as the
    values and GO id as the names}
  \item{goNEvi}{\code{goNEvi} a list of named vectors with GO ids as
    values for vectors and evidence code as names for vector values}
  \item{vector}{\code{vector} a vector that is going to be converted to
    a list using \code{as.list}}
  \item{vectNames}{\code{vectNames} a vector of character of string for
    the names of \code{vector} that is going to converted to a list}
  \item{srcObjs}{\code{srcobjs} a list that contains objects of the
    pubRepo class}
  \item{lazyLoad}{\code{lazyLoad} a boolean indicating whether a lazy
    load database will be created}
}
\details{
  These functions are the results of an effort to make data package
  building easier for urers. As the results, users may not have great
  power controlling the process or inputs. Additionally, some of the
  built in functions that figure out the urls for source data may fail
  when maintainers of the data source web sites change the name,
  structure, ect of the source data. When such event occurs, users may
  have to follow the instructions contained in a vignette named
  AnnBuilder to build data packages.

  \code{\link{getBaseParsers}} figures out which of the built in parsers
  to use to parse the source data based on the type of the mappings done
  for the probes.

  \code{\link{createEmptyDPkg}} creates an empty package with the
  required subdirectories for data to be stored.

  \code{\link{getMultiColNames}} figures out what data elements for
  annotation have many to one relations with a probe. The many parts are
  separated by a separator in parsed annotation data.

  \code{\link{getUniColNames}} figures out what data elements for
  annotation have one to one relations with a probe.

  \code{\link{getTypeColNames}} figures out what data elements for
  annotation have many to one relations with a probe and additional
  information appended to the end of each element following a
  separate. The many parts are also separated by a separator in parsed
  annotation data.

  \code{splitEntry} splits entries by a separator.
  
  \code{twoStepSplit} splits entries by the separator specified by sep
  and the descriptive information of each element by eleSep.
}
\value{
  \code{\link{getBaseParsers}} returns a named vector for the names of
  the parsers to use to parse the source data.

  \code{\link{getDirContent}} returns a vector of chracter strings for
  the content of a directory of interests.

  \code{\link{getMultiColNames}} returns a vector of character strings.
  
  \code{\link{getUniColNames}} returns a vector of character strings.
  
  \code{\link{getTypeColNames}} returns a vector of character strings.

  \code{splitEntry} returns a vector of character strings.

  \code{twoStepSplit} returns a named vector of character strings. The
  names are the descriptive information appended to each element by
  \code{eleSep} 
}
\references{ABPrimer and AnnBuilder vignettes}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{GOPkgBuilder}},\code{\link{KEGGPkgBuilder}}}
\examples{
# Create a temporary directory for the data
myDir <- tempdir()
# Create a temp base data file
geneNMap <- matrix(c("32468_f_at", "D90278", "32469_at", "L00693",
                   "32481_at", "AL031663", "33825_at", " X68733",
                   "35730_at", "X03350", "36512_at", "L32179",
                   "38912_at", "D90042", "38936_at", "M16652",
                   "39368_at", "AL031668"), ncol = 2, byrow = TRUE)
write.table(geneNMap, file = file.path(myDir, "geneNMap"),
sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# Urls for truncated versions of source data
mySrcUrls <- c(LL =
               "http://www.bioconductor.org/datafiles/wwwsources/Tll_tmpl.gz", UG = "http://www.bioconductor.org/datafiles/wwwsources/Ths.data.gz", 
GO = "http://www.bioconductor.org/datafiles/wwwsources/Tgo.xml")
# Create temp files for other sources
temp <- matrix(c("32468_f_at", NA, "32469_at", "2",
                   "32481_at", NA, "33825_at", " 9",
                   "35730_at", "1576", "36512_at", NA,
                   "38912_at", "10", "38936_at", NA,
                   "39368_at", NA), ncol = 2, byrow = TRUE)
write.table(temp, file = file.path(myDir, "srcone"), sep = "\t",
quote = FALSE, row.names = FALSE, col.names = FALSE)
temp <- matrix(c("32468_f_at", NA, "32469_at", NA,
                   "32481_at", "7051", "33825_at", NA,
                   "35730_at", NA, "36512_at", "1084",
                   "38912_at", NA, "38936_at", NA,
                   "39368_at", "89"), ncol = 2, byrow = TRUE)
write.table(temp, file = file.path(myDir, "srctwo"), sep = "\t",
quote = FALSE, row.names = FALSE, col.names = FALSE)
otherMapping <- c(srcone = file.path(myDir, "srcone"),
srctwo = file.path(myDir, "srctwo"))
# Runs only upon user's request
if(interactive()){
ABPkgBuilder(baseName = file.path(myDir, "geneNMap"),
srcUrls = mySrcUrls, baseMapType = "gb", otherSrc = otherMapping,
pkgName = "myPkg", pkgPath = myDir, organism = "Homo sapiens", version =
"1.1.0", makeXML = TRUE, author = c(author = "myname", maintainer =
"myname@myemail.com")) 
# Output files
list.files(myDir)
# Content of the data package
list.files(file.path(myDir, "myPkg"))
list.files(file.path(myDir, "myPkg", "data"))
list.files(file.path(myDir, "myPkg", "man"))
list.files(file.path(myDir, "myPkg", "R"))
unlink(file.path(myDir, "myPkg"), TRUE)
unlink(file.path(myDir, "myPkg.xml"))
unlink(file.path(myDir, "myPkgByNum.xml")) 
}
unlink(c(file.path(myDir, "geneNMap"), file.path(myDir, "srcone"),
file.path(myDir, "srctwo")))
}
\keyword{manip}

