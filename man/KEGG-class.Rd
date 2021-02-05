\name{KEGG-class}
\docType{class}
\alias{KEGG-class}
\alias{KEGG}
\alias{findIDNPath}
\alias{mapLL2ECNPName}
\alias{parseKEGGGenome}
\alias{KEGG,KEGG-method}
\alias{findIDNPath,KEGG-method}
\alias{mapLL2ECNPName,KEGG-method}
\title{Class "KEGG" a sub-class of pubRepo to get/process pathway and
  enzyme information}
\description{This class is a sub-class of pubRepo with source specific
  functions to get/process data from KEGG
  {\url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}} to obtain pathway and
  emzyme information for genes}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("KEGG", ...)}.
    A constructor (KEGG) is available and should be used to instantiate
    objects of this class 
}
\section{Slots}{
  \describe{
    \item{\code{organism}:}{Object of class \code{"character", from
	class "UG"} a character string for the organism of concern}
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"UG"} a character string for the url where source data are
      stored (\url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}) at the
      time of coding}
    \item{\code{parser}:}{Object of class \code{"character", from class
	"UG"} not in use}
    \item{\code{baseFile}:}{Object of class \code{"character", from
	class "UG"} not in use}
  }
}
\section{Extends}{
Class \code{"UG"}, directly.
Class \code{"pubRepo"}, by class "UG".
}
\section{Methods}{
  \describe{
    \item{findIDNPath}{\code{signature(object = "KEGG")}: Finds the
      mappings between KEGG ids and pathway names}
    \item{mapLL2ECNPName}{\code{signature(object = "KEGG")}: Maps
      LocusLink ids to enzyme ids and pathway names}
  }
}
\references{\url{www.genome.ad.jp/kegg/}}
\author{Jianhua Zhang}

\seealso{\code{\link{pubRepo-class}}, \code{\link{UG-class}} }

\examples{
\dontrun{
# The url (\url{ftp://ftp.genome.ad.jp/pub/kegg/pathways}) may change but
# was correct at the time of coding
url <-  getSrcUrl("KEGG")
kegg <- KEGG(srcUrl = url, organism = "human")
pathNEnzyme <- mapLL2ECNPName(kegg)
}}
\keyword{classes}
