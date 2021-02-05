\name{MeSHParser}
\alias{MeSHParser}
\alias{setVars}
\title{Function to parse the XML data file form MeSH }
\description{
  Given the name of a local version of the XML file from MeSH, this
  function parses the file and returns a list of environment objects
  containing the subtracted data.
}
\usage{
MeSHParser(mesh)
setVars()
}

\arguments{
  \item{mesh}{\code{mesh} a character string for the name of a local
    version of the XML data file avaiable for downloading from MeSH}
}
\details{
  Due to security reasons at the servers end, the source XML file has to
  be downloaded from MeSH and stored locally. MeSHParser reads the file
  to subtract data.
}
\value{
  This function returns a list of environment objects.
  \item{treenum}{a vector that contains mappings between Descriptor
    unique ids and their corresponding tree number assigned by MeSH}
  \item{scopenote}{a vector that contains mappings between Descriptor
    unique ids and their corresponding notes provided by MeSH}
  \item{qualifier}{a vector contains mappings between Descriptor and
    corresponding qualifier headings}
  \item{concept}{a vector contains mappings between the headings of Descriptor
    and corresponding Concepts belonging to the Descriptor}
  \item{term}{a vector contains mappings between the headings of
    Concepts and the corresponding Terms belonging to the Concepts}
  \item{heading}{a vector contains mappings between the unique MeSH ids
    and their corresponding headings}
}
\references{http://www.nlm.nih.gov/mesh/meshhome.html}
\author{Jianhua Zhang}

\keyword{manip}
