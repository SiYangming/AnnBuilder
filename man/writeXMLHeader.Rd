\name{writeXMLHeader}
\alias{writeXMLHeader}

\title{A function to write header information to an XML file.}
\description{
  This function writes to the Attr node of an annotate XML file. 
}
\usage{
writeXMLHeader(outName, fileCol, name, version, organism="human")
}

\arguments{
  \item{outName}{A character string for the name of the XML file to store
    the generated mata-data.}
  \item{fileCol}{A vector of character strings for the names of data
    columns in the original file that is going to be used to produce the
    Data node of the XML file.}
  \item{name}{A character string for an internal name that is normally
    the target of the annotation (e. g. U95 for the u95 chip).}
  \item{version}{A character string or number for the version of the
    system that produces the XML file.}
  \item{organism}{A character string for the name of the organism of
    interests} 
}
\details{
  The XML file produced has an Attr node to hold the header
  information. The Attr node contains a Target node for the internal
  name, a DataMade node to date the file when it is made, one to many
  SourceFile nodes for names of the source files used for annotation,
  and one to many Element nodes for names of the data items the Data
  node of the XML will contain.
}
\value{
  This function does not return any value.
}
\references{\url{http://www.bioconductor.org/datafiles/dtds/annotate.dtd}}
\author{Jianhua (John) Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R.}

\seealso{\code{\link{fileToXML}}}

\examples{
\dontrun{
makeSrcInfo()
#Write the header to a temp file
writeXMLHeader(outName = "try.xml", fileCol = c("AFFY", "LOCUSID",
"ACCNUM"), name = "Not Real", version = "0.5", organism = "human")
# View the header
readLines("try.xml")
# Clearn up
unlink("try.xml")
}
}
\keyword{manip}







