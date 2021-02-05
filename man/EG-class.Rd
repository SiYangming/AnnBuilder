\name{EG-class}
\docType{class}
\alias{EG-class}
\alias{parseData,EG-method}
\alias{EG}

\title{Class "EG" handles data provided by Entrez Gene}
\description{Entrez Gene contains data that were previously provided by
  LocusLink. The EG class represents objects that contains the needed
  information for gettting and processing the data}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EG", ...)}.
    A constructor (\code{\link{EG}} is available and should be used to
    instantiate objects of EG} 
}
\section{Slots}{
  \describe{
    \item{\code{accession}:}{Object of class \code{"character"} for the
      name of the file containing mappings between GenBank accession
      numbers and Gene ids}
    \item{\code{info}:}{Object of class \code{"character"} for the name
      of the file containing mappings between Gene ids and symbol,
      chromosome number for genes, cytoband information, and gene name}
    \item{\code{go}:}{Object of class \code{"character"} for the name of
    the file containing mappings between Gene ids and GO information}
    \item{\code{pubmed}:}{Object of class \code{"character"} for the
      name of the file containing mappings between Gene ids and PubMed ids}
    \item{\code{refseq}:}{Object of class \code{"character"} for the
      name of the file containing mappings between Gene ids and RefSeq ids}
    \item{\code{unigene}:}{Object of class \code{"character"} for the
      name of the file containing the mappings between Gene ids and
      UniGene ids}
    \item{\code{mim}:}{Object of class \code{"character"} for the name
      of the file containing mappings between Gene ids and OMIM ids}
    \item{\code{srcUrl}:}{Object of class \code{"character"} for the
      root URL where the aforementioned files reside}
    \item{\code{parser}:}{Object of class \code{"character"} for the
      name of a Perl parser that will be used to parse the source file}
    \item{\code{baseFile}:}{Object of class \code{"character"} for the
      name of the base file that contains mappings between probe ids and
    a public database ids that will be used to map probe ids to
    annotation data contained in a source file}
    \item{\code{built}:}{Object of class \code{"character"} for build
      information of the source file}
    \item{\code{fromWeb}:}{Object of class \code{"logical"} for
      inficating whether the source file should be accessed through the
      web or locally}
  }
}
\section{Extends}{
Class \code{"pubRepo"}, directly.
}
\section{Methods}{
  \describe{
    \item{parseData}{\code{signature(object = "EG")}: A method to parse
      a source file using a specified parser}
  }
}
\author{Jianhua Zhang}

\seealso{
  \code{\link{pubRepo-class}}
}
\keyword{classes}
