\name{athPkgBuilder}
\alias{athPkgBuilder}
\alias{getOneMap}
\alias{procPMIDData}
\alias{getSrcObjs4Ath}
\alias{readAthData}
\alias{mergeDupMatByFirstCol}

\title{Functions that build annotation packages for Arbidopsis}
\description{
  These functions are implemented specifically for building annotation
  data pckages for arabidopsis using the Arabidopsis information source
  (TAIR). 
}
\usage{
athPkgBuilder(baseName, pkgName, pkgPath,
                   fileExt = list(estAssign =
                   "Genes/est_mapping/est.Assignment.Locus",
                   seqGenes = "Genes/TAIR_sequenced_genes",
                   go = "Ontologies/Gene_Ontology/ATH_GO_GOSLIM.20050507.txt",
                   aliases = "Genes/gene_aliases.20041105",
                   pathway = "Pathways/aracyc_dump_20050412",
                   pmid = "Genes/Gene_Anatomy/ATH_Anatomy.20040209.txt"),
                   ncols = list(estAssign = 7, seqGenes = 4, go = 12,
                   aliases = 4, pathway = 4, pmid = 15),
    cols2Keep = list(estAssign = c(3, 6, 7), seqGenes = c(1, 3),
                   go = c(1, 5, 9), aliases = c(1, 2),
                   pathway = c(1, 3, 4), pmid = c(6, 11)),
    colNames = list(estAssign = c("CHRLOC", "ORI",  "ACCNUM"),
                   seqGenes = c("ACCNUM", "CHR"),
                   go = c("ACCNUM", "GO", "EVID"),
                   aliases = c("ACCNUM", "SYMBOL"),
                   pathway = c("PATH", "ENZYME", "ACCNUM"),
                   pmid = c("PMID", "ACCNUM")),
    version,
    author,
    baseUrl="ftp://ftp.arabidopsis.org/home/tair/", lazyLoad=TRUE)
getOneMap(map, keyCol)
procPMIDData(pmid)
getSrcObjs4Ath()
readAthData(baseUrl, ext, col2Keep, colNames, ncols)
mergeDupMatByFirstCol(dupMat, sep = ";")
}

\arguments{
  \item{baseName}{\code{baseName} a character string for the name of the
    base file to be used to build an annotation data package. The base
    file is assumed to have two columns with the first one being probe
    ids and second one being the corresponding TAIR locus ids.}
  \item{pkgName}{\code{pkgName} a character string for the name of the
    data package to be built}
  \item{pkgPath}{\code{pkgPath} a character string for the path to a
    directory where the data package to be built will be stored}
  \item{fileExt}{\code{fileExt} a list of character strings for the
    extension to be appended to a base url to form a complete url for a
    desired source data file stored at TAIR's ftp site. Some of the
    names given as default will change with time and need to be updated}
  \item{ncols}{\code{ncols} an integer indicating the total number of
    columns of a given source data file}
  \item{cols2Keep}{\code{cols2Keep} a vector of integers indicating
    which of the columns of a given source data file will be retained
    when the source file is read}
  \item{colNames}{\code{colNames} a vector of character strings for the
    names of the columns of the source file to be retained}
  \item{version}{\code{version} a character string for the version
    number of the data package to be built}
  \item{author}{\code{author} a list of character stirngs with an author
    and maintainer element for the name and email address of the author}
  \item{baseUrl}{\code{baseUrl} a character string for the base url to
    TAIRs ftp site, The default is
    \url{ftp://tairpub:tairpub@ftp.arabidopsis.org/home/tair/}}
  \item{map}{\code{map} a matrix containing mappings between probe ids
    and annotation data}
  \item{keyCol}{\code{keyCol} an integer or character string for the
    name of the column in a matrix that contains the keys based on which
    data in the other columns will be merged for duplicated keys}
  \item{pmid}{\code{pmid} a matrix containing mappings between probe ids
    and PubMed ids regarding genes represented by the probe ids}
  \item{ext}{\code{ext} a single string version of \code{fileExt}}
  \item{dupMat}{\code{dupMat} a matrix with duplicating values for
    entries in a column defined as keys}
  \item{sep}{\code{sep} a character string for separator to be used when
    values in a matrix are merged based on keys contained in another
    columns}
  \item{col2Keep}{\code{col2Keep} a vector of integers indicating which
    of the column of a data file will be kept when a file is read}
  \item{lazyLoad}{\code{lazyLoad} a boolean indicating whether a lazy
    load database will be created}

}
\details{
  The annotation data will be extracted from various sources that may
  change in both names and contents. The default values provided were
  correct at the time of implementation but may need updating when the
  function is actually used.
}
\value{
  The main function athPkgBuilder returns invisible()
}
\references{\url{http://www.arabidopsis.org}}
\author{Jianhua Zhang}

\seealso{\code{\link{ABPkgBuilder}}}
\examples{
# No example is provided due to the length of time required to build a package
}
\keyword{manip}


