\name{pfamDB}
\alias{getPfamDb}
\alias{closePfamDb}
\alias{setPfamDb}
\alias{getPfam}
\title{Database Utilities for PFAM package}
\description{
  Utilities to deal with the database in the PFAM package
}
\usage{
getPfamDb()
closePfamDb()
setPfamDb(con)
getPfam(pfamCon, FROM, TO, TABLE, FOR)
}
\arguments{
  \item{pfamCon}{The default connection to the RSQLite PFAM database}
  \item{con}{A new connection instead of pfamCon}
  \item{FROM}{The filed in the database which the mapping from}
  \item{TO}{The filed in the database which the mapping to}
  \item{TABLE}{The table in the database used to create the mapping}
  \item{FOR}{The selected value in the \code{FROM} field which the
    mapping focuses on}
}
\value{
  \code{getPfam} will return a data frame.  There is no return values
  for \code{getPfamDb}, \code{closePfamDb}, and \code{setPfamDb}.  
}
\details{
  \code{getPfamDb} will initialize the global connection to the PFAM
  database file, \code{closePfamDb} will close the global connection,
  and \code{setPfamDb} could set a new connection to the global
  connection.
  \code{getPfam} will send SQL quey to the database by the given
  arguments---\code{FROM}, \code{TO}, \code{TABLE}, and \code{FOR}---to
  retrieve the whole or the part of a table.  
}
\author{Ting-Yuan Liu}
\keyword{datasets}
