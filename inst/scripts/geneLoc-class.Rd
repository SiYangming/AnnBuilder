\name{geneLoc-class}
\docType{class}
\alias{chrom}
\alias{endLoc}
\alias{probeId}
\alias{print}
\alias{startLoc}
\alias{strand}
\alias{inExon}
\alias{geneLoc-class}
\alias{chrom,geneLoc-method}
\alias{endLoc,geneLoc-method}
\alias{probeId,geneLoc-method}
\alias{print,geneLoc-method}
\alias{startLoc,geneLoc-method}
\alias{strand,geneLoc-method}

\title{Class "geneLoc", a class of geneLocs to represent gene location data}
\description{A class of geneLocs is designed to represent data for a
  given length of sequence of on a chromosome. The class contains five
  slots that are of importance for a given length of chromosome sequence}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("geneLoc", probeId
= ..., chrom = ..., strand = ..., startLoc = ..., endLoc = ...)} 
}
\section{Slots}{
  \describe{
    \item{\code{probeId}:}{Object of class \code{"character"} a
      character string for the probe id corresponding the sequence}
    \item{\code{chrom}:}{Object of class \code{"character"} a character
      string for the chromosome number the sequence belongs to}
    \item{\code{strand}:}{Object of class \code{"character"} a character
      string for the strand of the chromosome that can either be - or +}
    \item{\code{startLoc}:}{Object of class \code{"numeric"} an integer
      for the starting location of the sequence on a chromosome}
    \item{\code{endLoc}:}{Object of class \code{"numeric"} an integer
      for the ending location of the sequence on a chromosome}
  }
}
\section{Methods}{
  \describe{
    \item{chrom}{\code{signature(object = "geneLoc")}: get method for
      chromosome number}
    \item{endLoc}{\code{signature(object = "geneLoc")}: get method for
      ending location of the sequence}
    \item{probeId}{\code{signature(object = "geneLoc")}: get method for
      probe id}
    \item{print}{\code{signature(x = "geneLoc")}: print method for the object}
    \item{startLoc}{\code{signature(object = "geneLoc")}: get method for
      stating location of the sequence}
    \item{strand}{\code{signature(object = "geneLoc")}: get method for strand}
  }
}
\references{\url{http://www.genome.ucsc.edu}}
\author{Jianhua Zhang}

\seealso{
  \code{\link{exonPkgBuilder}}
}
\examples{
require("AnnBuilder") || stop("Does not run without AnnBuilder")
gene <- new("geneLoc", probeId = "2", chrom = "12", strand = "-",
startLoc = 12345, endLoc = 12456)
probeId(gene)
startLoc(gene)
}
\keyword{classes}
