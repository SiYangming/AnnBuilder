\name{HG-class}
\docType{class}
\alias{HG-class}
\alias{HG}
\alias{readData,HG-method}

\title{Class "HG" a class to represent HomoloGene data source}
\description{Objects of HG contains the url, build information,
  ... about the HomoloGene data that will be used to build a homolgy
  data package}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("GO", ...)}.
    A constructor (\code{\link{HG}} is available and should be used to
    instantiate objects of HG).
}
\section{Slots}{
  \describe{
    \item{\code{srcUrl}:}{Object of class \code{"character"} a character
    string for the url or path of a source file to be used}
    \item{\code{parser}:}{Object of class \code{"character"} a character
    string for the name of a parser to be used to parse the source
    data. Not applicable to HG objects}
    \item{\code{baseFile}:}{Object of class \code{"character"} Not
      applicable to HG objects}
    \item{\code{built}:}{Object of class \code{"character"} a chracter
      string for the build information about the source file}
    \item{\code{fromWeb}:}{Object of class \code{"logical"} a boolean
      indicating whether \code{srcUrl} is a url to a source file or the
      path to a locally stored file}
  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
  \describe{
    \item{readData}{\code{signature(object = "HG")}: a function that
      reads the homoloGene data}
  }
}
\references{\url{http://www.ncbi.nlm. nih.gov/query?db = homology}}
\author{Jianhua Zhang}

\seealso{\code{\link{pubRepo-class}}
}
\keyword{classes}
