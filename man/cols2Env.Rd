\name{cols2Env}
\alias{cols2Env}
\alias{matchAll}
\alias{matchOneRow}

\title{Creates a environment object using data from two columns of a matrix}
\description{
  Given a matrix with two columns, this function creates an environment
  object with values in one of the specified columns as keys and those in the
  other column as values.
}
\usage{
cols2Env(cols, colNames, keyColName = colNames[1], sep)
matchAll(cols, keyColName)
matchOneRow(cols, keyColName, sep = ";")
}

\arguments{
  \item{cols}{\code{cols} a matrix with two columns}
  \item{colNames}{\code{colNames} a character string for the name of the
    column whose values will be used for the keys of the environment
    object to be created}
  \item{keyColName}{\code{keyColName} a character string for the name of
    the column whose values will be the corresponding values for keys of
    the environment object to be created}
  \item{sep}{\code{sep} a character for the separators used to separate
    entries that have multiple values}
}
\details{
  The matrix or matrix convertible object passed to cols2Env must have
  two coloumns with one intended to be used as the key and the other be the
  value.

  Cells in either or both columns may have multiple values separated by a
  separator (e.g. "a;b", "1;2;3") making the mapping between keys and
  the corresponding values not a straitforward operation. cols2Env gets
  all the unique values from the key column by spliting them and maps
  values to each of them.

  \code{\link{cols2Env}} calls \code{\link{matchAll}} that in turn calls
  \code{\link{matchOneRow}} to first split entries and then map entries
  in the two coloumns on one to one bases. Unique keys in the column
  defined as the key column will be assigned a vector containing all the
  values corresponding the keys in the environment to return. 
}
\value{
  This function returns an environment object with key and value pairs
}

\author{Jianhua Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{ABPkgBuilder}}}

\examples{
dataM <- matrix(c("a;b", "1;2;3", "a;b", "4;5", "c", "6;7", "b;a",
"6;7;8"), ncol = 2, byrow = TRUE)

temp <- AnnBuilder:::cols2Env(dataM, c("key", "value"), keyColName = "key")

dataM
mget(ls(temp), temp)
}
\keyword{manip}






