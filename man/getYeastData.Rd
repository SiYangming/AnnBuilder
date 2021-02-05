\name{getYeastData}
\alias{getYeastData}
\alias{readBadData}
\alias{findNumCol}
\title{Functions to get/process yeast genome data}
\description{
  These functions extract data from the yeast genome web site based on a
  set of arguments.
}
\usage{
getYeastData(url, extenName, cols2Keep, sep)
readBadData(url, sep)
findNumCol(fewLines, sep)
}
\arguments{
  \item{url}{\code{url} a character string for the url where yeast data
    are stored}
  \item{extenName}{\code{extenName} a character string for the name of
    the data file of interest. The name can be a file name or with
    subdirectory names under "url"}
  \item{cols2Keep}{\code{cols2Keep} a vector of index for the columns to
  be extracted from the data file}
  \item{sep}{\code{sep} a character string for the separater used to
    separate data columns in the data file}
  \item{fewLines}{\code{fewLines} a set of character strings separated
    by a new line that is going to be used to determine how many data
    columns each line has}
}
\details{
  The yeast genome web site has files stored in or in subdirectories of
  \url{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/} that can
  be downloaded. \code{\link{getYeastData}} extracts data from a given
  file. The functions are used by an object of \code{\link{YG-class}} to
  extract data.

  Some of the data in the web site may not be well fomatted (e.g. with
  missing columns). \code{\link{readBadData}} deals with these type of
  data files.

  \code{\link{findNumCol}} figures out how many data columns a file
  contains based on a few entries from that file.
  
}
\value{
  \code{\link{getYeastData}} returns a matrix containing data.

  \code{\link{readBadData}} returns a matrix.

  \code{\link{findNumCol}} returns an integer.
}
\references{\url{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/}}
\author{Jianhua Zhang}

\seealso{\code{\link{YG-class}}}
\examples{
\dontrun{
# Url may change but was correct at the time of coding
url <- "ftp://genome-ftp.stanford.edu/pub/yeast/data_download/"
temp <- getYeastData(url, "chromosomal_feature/SGD_features.tab",
                         cols2Keep = c(6, 1), sep = "\t")
}
}
\keyword{manip}

