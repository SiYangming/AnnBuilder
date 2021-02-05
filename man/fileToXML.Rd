\name{fileToXML}
\alias{fileToXML}

\title{A function to convert a text file to XML.}
\description{
  This function takes a text file and then converts the data contained
  by the file to an XML file. The XML file contains an Attr and a Data
  node. The Attr node contains mata-data and the Data node contains real
  data from the original file. 
}
\usage{
fileToXML(targetName, outName, inName, idColName, colNames,
multColNames, typeColNames,  multSep = ";", typeSep = ";", fileSep =
"\t", header = FALSE, isFile = TRUE, organism = "human",version = "1.0.0")
}

\arguments{
  \item{outName}{\code{outName}A character string for the name of xml file to be
    produced. If the name does not contain a full path, the current
    working directory will be the default}
  \item{inName}{\code{inName} A character string for the name of the input file to be
    written to an XML document}
  \item{idColName}{\code{idColName} A character string for the name of
    the column in the input file where ids of the target of annotation are}
  \item{colNames}{\code{colNames} A vector of character strings for the name of data
    columns in the original file.}
  \item{targetName}{\code{targetName} A character string that will be used as an internal name
    for the meta-data to show the target of the annotation (e.g. U95, U6800.}
  \item{version}{\code{version} A character string or number indicating the version of
    the system used to build the xml file.}
  \item{multColNames}{\code{multColNames} A vector of character strings for the name of data
    columns that may contain multiple items separated by a separator
    specified by parameter multSep. }
  \item{typeColNames}{\code{typeColNames} A vector of character strings for data columns in the
    original data that may contain type information append to the real
    data with a separator defined by parameter typeSep
    (e.g. "aGeneName;Officila").} 
  \item{multSep}{\code{mutlSep} A character string for the separator used to separate
    multiple data items within a data column of the original file.}
  \item{typeSep}{\code{typeSep} A character string for the separator used to separate
    the real data and type information within a column of the original data.}
  \item{fileSep}{\code{fileSep} A character string specifying how data columns are
    separated in the original file (e.g. sep = "\t" for tab delimited.}
  \item{organism}{\code{organism} A character string for the name of the organism of
    interests} 
  \item{header}{\code{header} A boolean that is set to TRUE if the original file has a
  header row or FALSE otherwise.}
  \item{isFile}{\code{isFile} A boolean that is set to TRUE if parameter fileName is
    the name of an existing file and FALSE if fileName is a R object
    contains the data}
}
\details{
  The original text file is assumed to have rows with columns separated
  by a separator defined by parameter sep. MultCol are used to define
  data columns that capture the one to many relationships between
  data. For example, a given AffyMetrix id may be associated with several
  GenBank accession numbers. In a data set with AffyMetrix ids as one of
  the data columns, the accession number column will be a element in
  multCol with a separator separating individual accession numbers
  (e.g. X00001,X00002,U0003... if the separator is a ",").

  As gene name and gene symbol can be "Official" or "Preferred", a type
  information is attached to a gene name or symbol that is going to be
  the value for attribute type in the resulting XML file
  (e.g. XXXX;Official if the separator is ";").
}
\value{
  This function does not return any value. The XML file will be stored
  as a file.
}
\references{\url{http://www.bioconductor.org/datafiles/dtds/annotate.dtd}}
\author{Jianhua (John) Zhang}

\seealso{\code{\link{ABPkgBuilder}}}

\examples{
# Create a text file
aFile <- as.data.frame(matrix(c(1:9), ncol = 3))

#Write to an XML file
if(interactive()){
    fileToXML("notReal", outName = "try.xml", inName = aFile, idColName =
    "AFFY", colNames = c("AFFY", "LOCUSID", "UNIGENE"), multColNames = NULL,
    typeColNames = NULL,  multSep = ";", isFile = FALSE)

    #Show the XML file
    readLines("try.xml")

    # Clearn up
    unlink("try.xml")
}
}
\keyword{manip}










