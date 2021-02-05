\name{cleanSrcObjs}
\alias{cleanSrcObjs}
\alias{getRepSourceNBuilt}
\alias{mapGO2Probe}
\alias{writeReverseMap}
\alias{writeAnnData2Pkg}
\alias{getAnnData}
\alias{getUniMappings}
\alias{getBaseFile}
\alias{getSrcObjs}
\alias{mapll2PathID}
\alias{mapLLNGB}
\alias{getLLNGBMap}
\alias{getRepList4Perl}
\alias{mapUGNGB}
\alias{getTaxid}

\title{Supporting function that may no of any other usese}
\description{
  Functions in this group are mainly for supporting purposes and may not
  be of any use outside the package they reside
}
\usage{  
getRepSourceNBuilt(name, object)
mapGO2Probe(eg, baseMapType)
writeReverseMap(annData, pkgName, pkgPath)
writeAnnData2Pkg(annData, pkgName, pkgPath)
getAnnData(srcObjs)
getUniMappings(baseName, eg, ug, otherSrc, baseMapType)
getBaseFile(baseName)
getSrcObjs(srcUrls, baseName, organism, baseMapType = c("gb", "ug",
"ll", "image", "refseq", "gbNRef"), fromWeb = TRUE)
cleanSrcObjs(srcObjs)
mapll2PathID(srcUrl, organism, exten = "gene_map.tab")
mapLLNGB(organism, pkgName, pkgPath, ugUrl = getSrcUrl("ug", organism),
egUrl = paste(getSrcUrl("eg"), "gene2accession.gz"), fromWeb = TRUE) 
getLLNGBMap(repList, what = "ll2gb")
mapUGNGB(organism, pkgName, pkgPath, ugUrl = getSrcUrl("ug", organism),
llUrl = getSrcUrl("ll"), fromWeb = TRUE)
getRepList4Perl(organism, ugUrl = getSrcUrl("ug", organism), llUrl =
getSrcUrl("ll"), fromWeb = TRUE) 
getTaxid(organism)
}

\arguments{
  \item{ugUrl}{\code{ugUrl} a character string for the url to the
    ftp site of UniGene}
  \item{llUrl}{\code{llUrl} a character string for the url to the ftp
    site of LocusLink}
  \item{egUrl}{\code{egUrl} a character string for the url to the ftp
    site of Entrez Gene}
  \item{repList}{\code{repList} a list with values to be used to replace
    contents in template files}
  \item{name}{\code{name} a character string for the name a data source
    to be used to build a data package}
  \item{object}{\code{object} an object that is a subclass of pubRepo}
  \item{ll}{\code{ll} an object of class LL}
  \item{ug}{\code{ug} an object of class ug}
  \item{baseMapType}{\code{baseMapType} a character string for the type
    of base map (e. g. gb, ug, ll, ...)}
  \item{annData}{\code{annData} a matrix derived from source data}
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be built}
  \item{pkgPath}{\code{pkgPath} a character string for the path where a
    new package will be built}
  \item{srcObjs}{\code{srcObjs} a list containing objects of class UG,
    LL, GO and so on}
  \item{baseName}{\code{baseName} a character string for the nam of a
    base file to be used to build a data package} 
  \item{otherSrc}{\code{otherSrc} a named vector for files contianing
    mappings between probe ids and LocusLink ids obtained by other
    sources}
  \item{srcUrls}{\code{srcUrls} a named vector for the urls to the
    source data to be used to build data packages}
  \item{srcUrl}{\code{srcUrl} a url for a source data file to be used}
  \item{fromWeb}{\code{fromWeb} a boolean indicating whether a source
    url is a real url or just the path to a locally stored file} 
  \item{organism}{\code{organism} a character string for the name of the
    organism of concern}
  \item{exten}{\code{exten} a character string for the extension to be
    appended to the end of a given url to make the url complete} 
  \item{what}{\code{what} a character string that can either be ll2gb or
    gb2ll}
  \item{eg}{\code{eg} an EG object}
}

\author{Jianhua Zhang}

\keyword{misc}



