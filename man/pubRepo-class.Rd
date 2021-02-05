\name{pubRepo-class} 
\docType{class}
\alias{pubRepo-class}
\alias{baseFile<-}
\alias{baseFile}
\alias{builtInfo}
\alias{downloadData}
\alias{parseData}
\alias{parser<-} 
\alias{parser}
\alias{readData}
\alias{srcUrl<-}
\alias{srcUrl}
\alias{baseFile<-,pubRepo-method}
\alias{baseFile,pubRepo-method}
\alias{downloadData,pubRepo-method}
\alias{parseData,pubRepo-method}
\alias{parser<-,pubRepo-method}
\alias{fromWeb<-,pubRepo-method}
\alias{parser,pubRepo-method}
\alias{readData,pubRepo-method}
\alias{srcUrl<-,pubRepo-method}
\alias{srcUrl,pubRepo-method}
\alias{builtInfo,pubRepo-method}
\alias{pubRepo}
\alias{fromWeb}
\alias{fromWeb<-}
\alias{fromWeb,pubRepo-method}

\title{Class "pubRepo" a generic class for downloading/parsing data
  provided by various public data repositories}
\description{This class provides basic functions to download/parse
  data from different public data repositories. More specific functions
  can be provided by extending this class to include source specific
  features}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pubRepo", ...)}.
    A constructor (\code{\link{pubRepo}} is provided and should be used
    to create objects of this class. 
}
\section{Slots}{
  \describe{
    \item{\code{srcUrl}:}{Object of class \code{"character"} a character
      string for the url of a data source from a public repository}
    \item{\code{parser}:}{Object of class \code{"character"} a character
      string for the name of a file that will be used as part of perl
      script to parse the source data. Parser is a segment of perl code
      containing instructions on how the source data will be processed
      and the content and format of the output}
    \item{\code{baseFile}:}{Object of class \code{"character"} a
      character string for the name of a file that will be used as the 
      base to process the source data. Data from the source that are
      related to elements in the base file will be extracted. baseFile
      is assumed to be a two column file with the first column being
      some type of arbitrary ids (e.g. Affymetrix probe ids) and the
      second column being the corresponding ids of a given public
      repository (e.g. GenBank accession numbers or UniGene ids)}
    \item{\code{built}:}{Object of class \code{"character"} a character
      string for the date or number a given source data were built}
    \item{\code{fromWeb}:}{Object of class \code{"boolean"} a boolean
      indicating whether the data will be read from a url or local file
      represented by \code{srcUrl}}
  }
}

\section{Methods}{
  \describe{
    \item{baseFile<-}{\code{signature(object = "pubRepo")}: Sets the
      value for baseFile}
    \item{baseFile}{\code{signature(object = "pubRepo")}: Gets the value
      for baseFile}
    \item{builtInfo}{\code{signature(object = "pubRepo")}: Gets the value
      for built}
    \item{downloadData}{\code{signature(object = "pubRepo")}: Downloads
      data from a data source defined by srcUrl}
    \item{parseData}{\code{signature(object = "pubRepo")}:
      DownLoads/parses data from a data source defined by srcUrl}
    \item{parser<-}{\code{signature(object = "pubRepo")}: Sets the value
      for parser}
    \item{parser}{\code{signature(object = "pubRepo")}: Gets the value
      for parser}
    \item{readData}{\code{signature(object = "pubRepo")}: Reads data
      using \code{\link{readLines}} from a data source defined by srcUrl}
    \item{srcUrl<-}{\code{signature(object = "pubRepo")}: Sets the value
      for srcUrl}
    \item{srcUrl}{\code{signature(object = "pubRepo")}: Gets the value
      for srcUrl}
    \item{fromWeb}{\code{signature(object = "pubRepo")}: Get the vlue
      for slot fromWeb}
    \item{fromWeb<-}{\code{signature(object = "pubRepo")}: Sets the value
      for slot fromWeb}
  }
}
\author{Jianhua Zhang}

\seealso{\code{\link{GO-class}}, \code{\link{KEGG-class}},
  \code{\link{LL-class}}, \code{\link{UG-class}}, \code{\link{GEO-class}}}

\examples{
\dontrun{
# Read a short test file from Bioconductor
test <- pubRepo(srcUrl =
"http://www.bioconductor.org/datafiles/wwwsources/TGene.txt", fromWeb = TRUE)
data <- readData(test)
}
}
\keyword{classes}


