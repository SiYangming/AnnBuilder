\name{getChroLocation}
\alias{getChroLocation}
\alias{getGPData}
\alias{gpLinkNGene}
\alias{getCytoList}
\alias{getCytoLoc}
\title{Functions to extract data from Golden Path}
\description{
  These functions are used by objects GP to extract chromosomal location
  and orientation data for genes using source files provided by Golden Path
}
\usage{
getChroLocation(srcUrl, exten = gpLinkNGene(), sep = "\t", fromWeb =
TRUE, raw = FALSE)
getGPData(srcUrl, sep = "\t")
gpLinkNGene(test = FALSE, fromWeb = TRUE)
getCytoList(data)
getCytoLoc(organism, srcUrl = paste(getSrcUrl("gp",organism), "/",
          "http://www.genome.ucsc.edu/goldenPath/mm4/database/" ))
}
\arguments{
  \item{srcUrl}{\code{srcUrl} a character string for the url where
    Golden Path source data are available}
  \item{exten}{\code{exten} a character string for the name of the file
    to be used for the extraction}
  \item{sep}{\code{sep} a character string for the separator used by the
    source file}
  \item{test}{\code{test} a boolean to indicate whether the process is
    in a testing mode}
  \item{fromWeb}{\code{fromWeb} a boolean to indicate whether the source
    data should be downloaded from the web or is a local file}
  \item{raw}{\code{raw} a boolean indicating whether chromosomal location
    data will be returned as a five column data frame with ID,
    Chromosome, strand, start, and end or a two column data with ID and
    processed chromosome location data}
  \item{organism}{\code{organism} a character string for the name of the
    organism of interest}
  \item{data}{\code{data} a data matrix}
}
\details{
  \code{\link{getChroLocation}} extracts chromosomal location data from
  a data file named refGene.

  \code{\link{getGPData}} Reads data from a source data file defined by
  srcUrl and returns them as a matrix.

  \code{\link{gpLinkNGene}} returns a correct link and gene data file
  names that will be used to get chromosomal location data.
}
\value{
  \code{\link{getChroLocation}} returns a matrix with five or two columns.

  \code{\link{getGPData}} returns a matrix.

  \code{\link{gpLinkNGene}} returns a named vector.
}
\references{\url{http://www.genome.ucsc.edu}}
\author{Jianhua Zhang}

\seealso{\code{\link{GP}}}
\examples{
\dontrun{
# Truncated versions of files stored in Bioconductor site are used
gpLinkNGene(test = FALSE)
temp <- getGPData(
"http://www.bioconductor.org/datafiles/wwwsources/Tlink.txt",
sep = "\t", ncol = 8, keep = c(3,7))
temp <- getChroLocation(
"http://www.bioconductor.org/datafiles/wwwsources/",
exten = gpLinkNGene(TRUE), sep = "\t")
}
}
\keyword{manip}





