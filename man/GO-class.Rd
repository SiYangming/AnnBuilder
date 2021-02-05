\name{GO-class}
\docType{class}
\alias{GO-class}
\alias{GO}
\alias{readData,GO-method}
\title{Class "GO" a class to handle data from Gene Ontology}
\description{This class is a sub-class of pubRepo that is implemented
  specifically to parse data from Gene Ontology. \code{\link{readData}}
  has been over written to process Gene Ontology data}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GO", ...)}.
    A constructor (\code{\link{GO}} is available and should be used to
    instantiate objects of GO} 
}
\section{Slots}{
  \describe{
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"pubRepo"} a character string for the url of the source data
      from Gene Ontology}
    \item{\code{parser}:}{Object of class \code{"character", from class
	"pubRepo"} not in use}
    \item{\code{baseFile}:}{Object of class \code{"character", from
	class "pubRepo"} not in use}
  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
  \describe{
    \item{readData}{\code{signature(object = "GO")}: Downloads/processes
      go\_xxx-termdb from Gene Ontology, where xxx is a date. If argument
      xml is set to be TRUE, the data file will be parsed and a matrix
      with three columns will be returned. The first column is for GO
      ids, second for the GO ids of its direct parents, and third for
      the ontology term defined by Gene Ontology. Otherwise, the data
      (not in xml form) will be read in using \code{\link{readLines}}}
  }
}
\references{\url{http://www.godatabase.org}}
\author{Jianhua Zhang}
\note{This class is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{pubRepo-class}}}
\keyword{classes}
