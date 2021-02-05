\name{GP-class}
\docType{class}
\alias{GP-class}
\alias{GP}
\alias{getStrand}
\alias{getStrand,GP-method}
\title{Class "GP" a sub-class of pubRepo to get/process data from GoldenPath}
\description{This class is a sub-class of pubRepo with source specific
  functions to get/process data from GoldenPath
  {\url{http://www.genome.ucsc.edu/goldenPath}} to obtain gene
  location and orientation data} 
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GP", ...)}.
    A constructor (GP) is available and should be used to instantiate
    objects of this class 
}
\section{Slots}{
  \describe{
    \item{\code{organism}:}{Object of class \code{"character", from
	class "UG"} s character string for the organism of concern}
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"UG"} a character string for the url where the source data
      are. As multiple data sources will be used, srcUlr in this case is
      the location where the source data are
      (e.g. \url{http://www.genome.ucsc.edu/goldenPath/14nov2002/database/})}  
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
    \item{getStrand}{\code{signature(object = "GP")}: Processes the
      refLink and refGene data files and returns a matrix with gene
      location and orientation data}
  }
}
\references{\url{http://www.genome.ucsc.edu}}
\author{Jianhua Zhang}

\seealso{\code{\link{pubRepo-class}}}

\examples{
# The example may take a few second to finish
\dontrun{
## The url (\url{ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/})
## was correct at the time of coding. Replace with a correct one if it
## is invalid 
url <- getSrcUrl("GP", organism = "human")
gp <- GP(srcUrl = url, organism = "human")
strand <- getStrand(gp)
}
}
\keyword{classes}
