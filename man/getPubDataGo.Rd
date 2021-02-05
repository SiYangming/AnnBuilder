\name{getPubDataGo}
\alias{getPubDataHomoloGene}
\alias{getPubDataLocusLink}
\alias{getPubDataUniGene}
\alias{getPubDataEntrezGene}
\alias{getPubDataGoldenPath}
\alias{getPubDataGo}
\alias{getPubDataYeastGenome}
\alias{getPubDataKegg}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to download public domain annotation data sources}
\description{
  These functions are intended to help create a local mirror of public
  domain annotation data sources.  They all depend on having wget
  available.
  
}
\usage{
getPubDataHomoloGene(baseUrl, passive)
getPubDataLocusLink(baseUrl, passive)
getPubDataUniGene(baseUrl, passive)
getPubDataEntrezGene(baseUrl, passive)
getPubDataGoldenPath(baseUrl, passive)
getPubDataGo(baseUrl, passive)
getPubDataYeastGenome(baseUrl, passive)
getPubDataKegg(baseUrl, passive)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseUrl}{The URL.  Note that for some sources this is a
  directory and for others it fully specifies a file we want to have
  available.}
  \item{passive}{logical.  If \code{TRUE}, pass the \code{-passive-ftp}
  flag to \code{wget}}
}
\author{Seth Falcon <sfalcon@fhcrc.org>}
\keyword{data}

