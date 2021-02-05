\name{YG-class}
\docType{class}
\alias{YG-class}
\alias{YG}
\alias{readData,YG-method}
\title{Class "YG" a sub-class of pubRepo that reads/downloads data from
  yeast genomic}
\description{This class is a sub-class ob pubRepo that has source
  specific functions to extract data from Yeast Genome ftp site
  (\url{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/})}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("YG", ...)}.
    A constructor (YG) is available and should be used to instantiate
    objects of this class
}
\section{Slots}{
  \describe{
    \item{\code{srcUrl}:}{Object of class \code{"character", from class
	"pubRepo"} a character string for the url where surce data are
      available
      (\url{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/} at
      the time of coding)}
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
    \item{readData}{\code{signature(object = "YG")}: Reads source data
      defined by argument extenName from the ftp site}
  }
}
\references{\url{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/}}
\author{Jianhua Zhang}

\seealso{\code{\link{pubRepo-class}}}

\examples{
\dontrun{
# Url may change but was correct at the time of coding
url <- "ftp://genome-ftp.stanford.edu/pub/yeast/data_download/"
# Creat a YG object
ygeno <- YG(srcUrl = url)
# Read the file named "chromosomal_feature.tab". Takes a few
# seconds to finish
data <- readData(ygeno,
                 "chromosomal_feature/chromosomal_feature.tab",
                 cols2Keep =c(6, 1), sep = "\t")
}
}
\keyword{classes}
