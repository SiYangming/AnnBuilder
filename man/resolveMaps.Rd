\name{resolveMaps}
\alias{resolveMaps}
\alias{getVote}
\alias{hasDelimit}
\alias{getUnified}
\alias{getNoDup}
\title{Functions to obtain unified mappings for a given set of ids using
various sources} 
\description{
  These functions are used to obtain unified mappings between
  two sets of ids based on the mappings available from different
  sources. Each source provide mappings between two sets of ids.
}
\usage{
resolveMaps(maps, trusted, srcs, colNames = NULL, outName = "", asFile = TRUE)
getVote(voters, sep = ";")
getUnified(voters)
getNoDup(voters)
hasDelimit(entry, deli = ";") 
}
\arguments{
  \item{maps}{\code{maps} a matrix with mappings for a set of key ids to
    another set of ids provided by different sources. The first column
    is assumed to be the key ids and the rest are mappings to another set
    of ids provided by different sources}
  \item{trusted}{\code{trusted} a vector of characters to indicate the
    column number of "maps" whose mappings are more reliable and should
    be used when there are conflicts among sources}
  \item{srcs}{\code{srcs} a vector of character strings for the names of
    columns that contain mappings from different sources}
  \item{colNames}{\code{colNames} a vector of character strings for the
    names of columns in "maps"}
  \item{outName}{\code{outName} a character string for the name of the
    file to contain the unified mappings}
  \item{asFile}{\code{asFile} a boolean to indicate whether the unified
    mappings will be saved as a file}
  \item{voters}{\code{voters} a vector containing mappings from
    different sources}
  \item{entry}{\code{entry} a character string to be checked for the
    existence of a separator}
  \item{deli}{\code{deli} a character string for a separator}
  \item{sep}{\code{sep} same as deli}
}
\details{
  Each source may have different mappings from the key ids to another
  set of ids. \code{\link{resolveMaps}} resolves the conflicts and
  derives a set of unified mappings based on the mappings provided from
  several sources.

  \code{\link{getVote}} resolves the mappings for a given key id and
  returns a vector with unified mapping and the number of sources that
  agree with the unified mapping.

  \code{\link{getUnified}} finds agreement among values in a vector
  passed. If some values agree, get the one agreed by most sources.

  \code{\link{getNoDup}} gets a value based on predefined rules when
  values from different sources do not agree.
  
  \code{\link{hasDelimit}} checks to see if a delimiter exists}
}
\value{
  \code{\link{resolveMaps}} returns a matrix with the first column being
  the key id set, second being the unified mappings to another id set,
  and third the total number of agreements found among sources.

  \code{\link{getVote}} returns a two element vector.

  \code{\link{getUnified}} returns a character string.

  \code{\link{getNoDup}} returns a character string.

  \code{\link{hasDelimit}} returns TRUE or FALSE.
}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{LL-class}}, \code{\link{UG-class}}}
\examples{
\dontrun{
maps <- matrix(c("id1", "a", "a", "b", "id2", "c","d", "c",
"id3", "e","e", "e", "id4", NA, "f", NA, "id5", "g", NA, "h", "id6", NA,
"NA", "i", "id7", NA, NA, NA), ncol = 4, byrow = TRUE)
unified <- resolveMaps(maps, c("srcll", "srcug"),
c("srcll", "srcug", "srcgeo"),
colNames = c("key1", "srcll", "srcug", "srcgeo"), outName = "",
asFile = FALSE)
}
}
\keyword{manip}



