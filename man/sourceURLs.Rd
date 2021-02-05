\name{sourceURLs}
\alias{sourceURLs}
\title{A data file contains urls for data available from various public
  repositories} 
\description{
  This data file is used by various objects (through
  \code{\link{getSrcUrl}}) to get the correct urls for various data
  sources to be processed. 
}
\details{
  sourceURLs[[XX]] will get the url for data source XX, where XX is a
  short name for a particular public data repository. Valid names
  include "LL" - LocusLink, "UG" - UniGene, "GP" - GoldenPath, "GO" -
  Gene Ontology, "KEGG" - Kyoto Encyclopedia of Genes and Genomes, "GEO"
  - Gene Expression Omnibus, and "YG" - Yeast Genome.
}
\author{Jianhua Zhang}
\note{The data file part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{pubRepo-class}}}
\examples{
data("sourceURLs", package="AnnBuilder")
sourceURLs[["KEGG"]]
}
\keyword{file}

