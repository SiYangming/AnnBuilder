\name{GOXMLParser}
\alias{GOXMLParser}
\title{Parse the Gene Ontology OBO XML data file}
\description{
  Generate R \code{environment} objects containing data parsed from the
  Gene Ontology (GO) XML data file.  The GO data file is available from
  \url{http://www.geneontology.org}.  This parser is designed to parse
  the \code{go_YYYYMM-termdb.obo-xml} file.
}
\usage{
GOXMLParser(fileName)
}
\arguments{
  \item{fileName}{Name of the XML file containing the GO source data in
  obo-xml format.}
}
\value{
  A list of \code{environment} objects representing the GO data
  structures.  The environments returned are:
  \item{TERM}{See GOTERM environment in the GO package}
  \item{BPPARENTS}{See the GOBPPARENTS environment in the GO data packagey}
  \item{MFPARENTS}{See the GOMFPARENTS environment in the GO data package}
  \item{CCPARENTS}{See the GOCCPARENTS environment in the GO data package}
  \item{BPCHILDREN}{See the GOBPCHILDREN environment in the GO data package}
  \item{MFCHILDREN}{See the GOMFCHILDREN environment in the GO data package}
  \item{CCCHILDREN}{See the GOCCCHILDREN environment in the GO data package}
  \item{OBSOLETE}{See the GOOBSOLETE environment in the GO data package}
  \item{BPOFFSPRING}{See the GOBPOFFSPRING environment in the GO data package}
  \item{MFOFFSPRING}{See the GOMFOFFSPRING environment in the GO data package}
  \item{CCOFFSPRING}{See the GOCCOFFSPRING environment in the GO data package}
  \item{BPANCESTOR}{See the GOBPANCESTOR environment in the GO data package}
  \item{MFANCESTOR}{See the GOMFANCESTOR environment in the GO data package}
  \item{CCANCESTOR}{See the GOCCANCESTOR environment in the GO data package}
}
\author{Chenwei Lin, John Zhang, Seth Falcon}
\keyword{manip}

