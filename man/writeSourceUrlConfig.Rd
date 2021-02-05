\name{writeSourceUrlConfig}
\alias{writeSourceUrlConfig}
\title{Create a source URL config file}
\description{
  After creating a local mirror of the public data sources, use this
  function to create a config file suitable for reading back into R
  using \code{readSourceUrlConfig}.
}
\usage{
writeSourceUrlConfig(file)
}
\arguments{
  \item{file}{where to write the config file}
}
%- maybe also 'usage' for other objects documented here.
\details{
  The KEGG URL is handled as a special case.
}
\author{S. Falcon}
\keyword{data}
