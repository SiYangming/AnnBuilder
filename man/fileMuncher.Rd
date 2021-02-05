\name{fileMuncher}
\alias{fileMuncher}
\alias{mergeRowByKey}
\title{Dynamically create a Perl script to parse a source file base on
  user specifications}
\description{
  This function takes a base file, a source file, and a segment of Perl
  script specifying how the source file will be parsed and the generates
  a fully executable Perl script that is going to be called to parse the
  source file. 
}
\usage{
fileMuncher(outName, baseFile, dataFile, parser, isDir = FALSE)
mergeRowByKey(mergeMe, keyCol = 1, sep = ";")
}

\arguments{
  \item{outName}{\code{outName} a character string for the name of the file
    where the parsed data will be stored}
  \item{baseFile}{\code{baseFile} a character string for the name of the
    file that is going to be used as the base to process the source
    file. Only data that are corresponding to the ids defined in the
    base file will be processed and mapped}
  \item{dataFile}{\code{dataFile} a character string for the name of the
    source data file}
  \item{parser}{\code{parser} a character string for the name of the
    file containing a segment of a Perl script for parsing the
    source file. An output connection to OUT that is for storing parsed
    data, an input connection to BASE for importing base file, and an
    input connection to DATA for reading the source data file are
    assumed to be open. parser should define how BASE and DATA will be
    used to extract data and then store them in OUT}
  \item{isDir}{\code{isDir} a boolean indicating whether dataFile is a
    name of a directory (TRUE) or not (FALSE)}
  \item{mergeMe}{\code{mergeMe} a data matrix that is going to be
    processed to merge rows with duplicating keys}
  \item{keyCol}{\code{keyCol} an integer for the index of the column
    containing keys based on which entries will be merged}
  \item{sep}{\code{sep} a character string for the separator used to
    separate multiple values}
}
\details{
  The system is assumed to be able to run Perl. Perl scripts generated
  dynamically will also be removed after execution.

  \code{\link{mergeRowByKey}} merges data based on common keys. Keys
  multiple values for a given key will be separated by "sep".
}
\value{
  \code{\link{fileMuncher}} returns a character string for the name of
  the output file

  \code{\link{mergeRowByKey}} returns a matrix with merged data.
}

\author{Jianhua Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{resolveMaps}}}

\examples{
if(interactive()){
path <- file.path(.path.package("AnnBuilder"), "scripts")
temp <- matrix(c("32469_f_at", "D90278", "32469_at", "L00693", "33825_at",
"X68733", "35730_at", "X03350", "38912_at", "D90042", "38936_at",
"M16652"), ncol = 2, byrow = TRUE)
write.table(temp, "tempBase", sep = "\t", quote = FALSE,
row.names = FALSE, col.names = FALSE)
# Parse a truncated version of LL\_tmpl.gz from Bioconductor
srcFile <-
loadFromUrl("http://www.bioconductor.org/datafiles/wwwsources/Tll_tmpl.gz")  
fileMuncher(outName = "temp", baseFile = "tempBase", dataFile = srcFile,
parser =  file.path(path, "gbLLParser"), isDir = FALSE)
# Show the parsed data
read.table(file = "temp", sep = "\t", header = FALSE)
unlink("tempBase")
unlink("temp")
}
}
\keyword{manip}


