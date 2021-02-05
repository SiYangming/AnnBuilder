\name{getSrcBuilt}
\alias{getSrcBuilt}
\alias{getLLBuilt}
\alias{getUGBuilt}
\alias{getUCSCBuilt}
\alias{getGOBuilt}
\alias{getKEGGBuilt}
\alias{getYGBuilt}
\alias{getHGBuilt}
\alias{getRefBuilt4HS}
\alias{getRefSeqBuilt}
\alias{getEGBuilt}
\title{Functions that get the built date or number of the source data
  used for annotation}
\description{
  Given a data source name and organism, the built date or number of the
  annotation source data will be returned. The built date or number is
  provided by the data source through its web site. 
}
\usage{
getSrcBuilt(src = "LL", organism = "Homo sapiens")
getLLBuilt(url = "http://www.ncbi.nlm.nih.gov/LocusLink/statistics.html")
getUGBuilt(organism, url = "ftp://ftp.ncbi.nih.gov/repository/UniGene")
getUCSCBuilt(organism)
getGOBuilt(url = "http://www.godatabase.org/dev/database/archive/latest")
getKEGGBuilt(url = "http://www.genome.ad.jp/kegg/kegg2.html")
getYGBuilt()
getHGBuilt()
getRefSeqBuilt(organism)
getRefSeqBuilt()
getEGBuilt()
}

\arguments{
  \item{src}{A character string for name of the data source. See details
    for valid names} 
  \item{organism}{A character string for the name of the organism of
    interests. See details for valid names}
  \item{url}{A character string for the url from which built information
    can be obtained}
}
\details{
  \code{getLLBuilt} finds the built data for LocusLink from the statistics
  page.
  
  \code{getUGBuilt} finds the built data for UniGene from the Xx.info file,
  where Xx is the short organism name (e.g. Hs for human)
  
  \code{getUCSCBuilt} finds the built data for the Human Genome Project from
  the folder for the latest release.
  
  \code{getGOBuilt} finds the built data for Gene Ontology from the timestamp
  for the -ont.xml.gz file.
  
  \code{getKEGGBuilt} finds the built data for KEGG from kegg2.html page
  (Release version and date)

  \code{YGBuilt} gets built information for Yeast Genome data.

  Valid data source names include LL - LocusLink, UG - UniGene, UCSC -
  the Human Genome Project, GO - Gene Ontology, KEGG - KEGG, YG - Yeast
  Genome.

  Valid organism name include human, mouse, rat, and yeast at this time.
}
\value{
  All functions return a string for the built information
}
\references{\url{http://www.ncbi.nlm.nih.gov/LocusLink/statistics.html},
  \url{ftp://ftp.ncbi.nih.gov/repository/UniGene},
  \url{http://www.godatabase.org/dev/database/archive/latest},
  \url{http://www.genome.ad.jp/kegg/kegg2.html},
  \url{ftp://ftp.ncbi.nih.gov/refseq/LocusLink/},
  \url{http://www.yeastgenome.org}}
\author{Jianhua Zhang}
\note{The functions are part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{getSrcUrl}}}

\examples{
\dontrun{
# Get built information for LocusLink
ll <- getSrcBuilt(src = "LL")
ug <- getSrcBuilt(src = "UG", organism = "Homo sapiens")
yg <- getYGBuilt()
ll
ug
yg
}}
\keyword{datasets}

