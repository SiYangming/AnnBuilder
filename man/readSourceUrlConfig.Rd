\name{readSourceUrlConfig}
\alias{readSourceUrlConfig}
\title{Read a data source URL config file}
\description{
  Read a data source URL config file, a simple text file with two named
  columns, name and url.  
}
\usage{
readSourceUrlConfig(file, urlPrefix)
}
\arguments{
  \item{file}{path containing names and URLs}
  \item{urlPrefix}{If present, this will be prepended to all URLs parsed
    in \code{file}}
}
\value{
  A named list of URLs.
}
\author{S. Falcon}
\keyword{data}

