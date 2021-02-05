\name{makeLLDB}
\alias{makeLLDB}
\title{Create a lazy loading database for package data files}
\description{
  This function processes the \code{*.rda} files in a package's data
  subdirectory and replaces them with a lazy load database.
}
\usage{
makeLLDB(packageDir, compress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{packageDir}{Path to the package source directory}
  \item{compress}{ If \code{TRUE}, compress the resulting lazy database.}
}
\details{
  The purpose is to create a lazy load database before INSTALL time.
  This makes installation of source packages much faster because the
  lazy database has been precomputed.

  We needed this because we want the meta data packages to have lazy
  load symantics for the data objects.  Users should be able to load a
  data package using \code{require} and then ask for any of the data
  environments by name.  We want lazy loading of these data sets because
  they tend to contain large environments which would take a long time
  to load if we did it at attach time.

  }
\value{
  This function is called for its side-effect: creating a lazy loading
  database for a package's data files.

  Note that this function is destructive in that it removed the data
  files (the \code{*.rda} files) after creating the lazy database.
}
\author{R. Gentleman}
\keyword{manip}

