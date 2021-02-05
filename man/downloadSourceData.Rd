\name{downloadSourceData}
\alias{downloadSourceData}
\title{Create a local mirror of annotation data sources}
\description{
  Uses wget to mirror relevant portions of publicly available annotation
  data sources.  The goal is to create a local mirror that can be served
  on your LAN to reduce network load when building multiple annotation
  data packaages.
}
\usage{
downloadSourceData(passive=FALSE)
}
\arguments{
  \item{passive}{logical.  If \code{TRUE}, pass the \code{-passive-ftp}
    flag to \code{wget}}
}
%- maybe also 'usage' for other objects documented here.
\details{
  The data files will be downloaded to the current working directory.
  The KEGG pathway data is a special case.  We download the current
  tarball of the pathway data, but it needs to be unpacked in
  \code{kegg/pathways}.

  On unix-like systems, the KEGG data will be unpacked automagically.
}
\author{S. Falcon}
\keyword{data}

