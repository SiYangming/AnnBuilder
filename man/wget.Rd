\name{wget}
\alias{wget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrapper for system wget}
\description{
  A convenience wrapper to download/mirror websites.  Relies upon 
  wget being available in PATH as it is called via \code{system}.
}
\usage{
wget(url, levels, accepts, passive=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url}{ The URL to get}
  \item{levels}{ Recursion depth, see wget man page and the
    \code{--level} option}
  \item{accepts}{ character vector.  Gets passed to wget as the value of
    the \code{--accept} option}
  \item{passive}{logical.  If \code{TRUE}, pass the \code{-passive-ftp}
    flag to \code{wget}}
}
\author{S. Falcon}
\keyword{data}
