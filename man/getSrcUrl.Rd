\name{getSrcUrl}
\alias{getSrcUrl}
\alias{getAllUrl}
\alias{getLLUrl}
\alias{getUGUrl}
\alias{getUCSCUrl}
\alias{getGOUrl}
\alias{getKEGGUrl}
\alias{readURL}
\alias{getGEOUrl}
\alias{getYGUrl}
\alias{getHGUrl}
\alias{getRefSeqUrl}
\alias{getEGUrl}
\title{Functions that find the correct url for downloading annotation data}
\description{
  Given a source data name and organism name, the url from which the
  source annotation data can be downloaded will be returned. 
}
\usage{
getSrcUrl(src, organism = "Homo sapiens", xml = TRUE, dateOnly = FALSE)
getAllUrl(organism)
getLLUrl()
getUCSCUrl(organism, downloadSite)
getUGUrl(organism)
getGOUrl(xml = TRUE, dateOnly = FALSE) 
getKEGGUrl()
readURL(url)
getGEOUrl()
getYGUrl()
getHGUrl()
getRefSeqUrl(organism)
getEGUrl()
}

\arguments{
  \item{src}{A character string for the name of the data source. See
    details for valid names}
  \item{organism}{A character string for the name of the organism of
    interests}
  \item{url}{A character string for the url where the source data can be
    downloaded}
  \item{dateOnly}{A boolean that is set to TRUE if only the built date
    of the data souce will be returned or TRUE if the source url will be
    returned}
  \item{xml}{A boolean indicating whether the XML format data file will
    be downloaded/processed}
  \item{downloadSite}{\code{downLoadSite} a character string for the url
    to the general downloading site for the human, mouse, and rat data}
}
\details{
  \code{getAllUrl} finds the urls for all the data source including
  LocusLink, UinGene, the Human Geneome Project, Gene Ontology, and
  KEGG.
  
  \code{getLLUrl} finds the url fro LocusLink.
  
  \code{getUCSCUrl} finds the url for the Human Genome Project.
  
  \code{getUGUrl} finds the url for UniGene.
  
  \code{getGOUrl} finds the url for Gene Ontology.
  
  \code{getKEGGUrl} finds the url for KEGG.

  \code{getGEOUrl} finds the url for GOE (the CGI script)

  \code{getYGUrl} gets the url to the ftp site where Yeast Genome data
  can be downloaded.

  Valid data source names include LL - LocusLink, UG - UniGene, UCSC -
  the Human Genome Project, GO - Gene Ontology, KEGG - KEGG, and YG -
  Yeast Genome.

  Valid organism name include human, mouse, rat, and yeast at this time.
}
\value{
  \code{getAllUrl} returns a vector of character strings and all the others
  return a character string for the url
}
\references{\url{"http://www.ncbi.nlm.nih.gov/LocusLink/statistics.html"},
  \url{"ftp://ftp.ncbi.nih.gov/repository/UniGene"},
  \url{"http://www.godatabase.org/dev/database/archive/latest"},
  \url{"http://www.genome.ad.jp/kegg/kegg2.html"},
  \url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/},
  \url{http://www.yeastgenome.org}}
\author{Jianhau Zhang}
\note{The functions are part of the Buiconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{getSrcBuilt}}}

\keyword{datasets}












