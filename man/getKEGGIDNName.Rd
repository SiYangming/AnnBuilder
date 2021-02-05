\name{getKEGGIDNName}
\alias{getKEGGIDNName}
\alias{getKEGGOrgName}
\alias{getLLPathMap}
\alias{mapll2EC}
\alias{parseEC}
\title{Functions to get/process pathway and enzyme data from KEGG}
\description{
  These functions extract pathway and enzyme data from KEGG
  \url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}. The functions are used
  by \code{\link{KEGG-class}}.
}
\usage{
getKEGGIDNName(object, exten = "/map\_title.tab")
getKEGGOrgName(name)
getLLPathMap(srcUrl, idNName, organism, fromWeb = TRUE)
mapll2EC(id, srcUrl, organism, fromWeb, sep = "\t")
parseEC(llNEC)
}
\arguments{
  \item{srcUrl}{\code{srcUrl} a character string for the url where
    source data are available}
  \item{object}{\code{object} a KEGG object with the slots filled with
    correct values}
  \item{exten}{\code{exten} a character string for data file name as an
    extension}
  \item{name}{\code{name} a character string for the name of the
    organism of concern. "human", "mouse", and "rat" are the valid
    values for now}
  \item{organism}{\code{organism} same as name}
  \item{idNName}{\code{idNName} a named vector normally obtained obtained
    by using function \code{\link{getKEGGIDNName}}}
  \item{sep}{\code{sep} a character string for the separators used to
    separator entries in a file}
  \item{llNEC}{\code{llNEC} a line of tab separated character strings
    with the first character string being a LocusLink id and second
    being the mapping enzyme (EC) names}
  \item{id}{\code{id} a character string for the KEGG id used for
    different pathway files}
  \item{fromWeb}{\code{fromWeb} a boolean to indicate whether a source
    data file will be read from a web site or locally}
}
\details{
  \code{\link{getKEGGIDNName}} read the data file "map\_title.tab" from
  KEGG to obtain the mappings between KEGG ids and pathway
  names.

  \code{\link{getKEGGOrgName}} takes the name for an organism and
  returns a short verion of the name used by KEGG for that organism.

  \code{\link{getLLPathMap}} maps LocusLink ids to pathway and enzyme
  names for an organism using various data files from KEGG.

  \code{\link{mapll2EC}} maps LocusLink ids to enzyme (EC) names for a
  given pathway.

  \code{\link{parseEC}} extracts enzyme data from a line of tab
  separated character strings to map a LocusLink id to enzyme (EC) names.
}
\value{
  \code{\link{getKEGGIDNName}} returns a named vector with KEGG ids
  being the names and pathway names being values.

  \code{\link{getKEGGOrgName}} returns a character string.

  \code{\link{getLLPathMap}} returns a list of two elements named "llec"
  and "llpathname". Each element is a matrix with mappings between
  LocusLink ids to enzyme or pathway names.

  \code{\link{mapll2EC}} returns a matrix with the first column being
  LocusLink ids and second enzyme (EC) names.

  \code{\link{parseEC}} returns two elements vector with the first
  element being a LocusLink id and second being the mapping enzyme (EC)
  names. 
}
\references{\url{www.genome.ad.jp/kegg/}}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{KEGG-class}}}
\examples{
\dontrun{
getKEGGOrgName("Homo sapiens")
# This group of code needs a while to finish
# Url may change but was correct at the time of coding
idNPath <- getKEGGIDNName(KEGG(organism = "Homo sapiens"))
temp <- getLLPathMap("ftp://ftp.genome.ad.jp/pub/kegg/pathways",
idNPath, "Homo sapiens")
temp <- mapll2EC("00010", "ftp://ftp.genome.ad.jp/pub/kegg/pathways",
"Homo sapiens", sep = "\t")
}
}
\keyword{manip}
