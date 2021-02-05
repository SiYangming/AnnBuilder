\name{yeastAnn}
\alias{yeastAnn}
\alias{getProbe2SGD}
\alias{procYeastGeno}
\alias{formatGO}
\alias{formatChrLoc}
\alias{getGEOYeast}
\alias{getYGExons}
\title{Functions to annotate yeast genom data}
\description{
  Given a GEO accession number for a yease data set and the extensions
  for annotation data files names that are available from Yeast Genom
  web site, the functions generates a data package with containing
  annoatation data for yeast genes in the GEO data set.
}
\usage{
yeastAnn(base = "", yGenoUrl,
                 yGenoNames =
                 c("literature_curation/gene_literature.tab",
                 "chromosomal_feature/SGD_features.tab",
                 "literature_curation/gene_association.sgd.gz"), toKeep =
                 list(c(6, 1), c(1, 5, 9, 10, 12, 16, 6), c(2, 5, 7)),
                 colNames = list(c("sgdid", "pmid"), c("sgdid",
                 "genename", "chr", "chrloc", "chrori", "description",
                 "alias"), c("sgdid", "go")), seps = c("\t", "\t",
                 "\t"), by = "sgdid")
getProbe2SGD(probe2ORF = "", yGenoUrl,
             fileName = "literature_curation/orf_geneontology.tab",
             toKeep = c(1, 7), colNames = c("orf", "sgdid"), sep = "\t",
             by = "orf")
procYeastGeno(baseURL, fileName, toKeep, colNames, seps = "\t")
getGEOYeast(GEOAccNum, GEOUrl, geoCols = c(1, 8), yGenoUrl) 
formatGO(gos, evis)
formatChrLoc(chr, chrloc, chrori)
getYGExons(srcUrl,
           yGenoName = "chromosomal_feature/intron_exon.tab", sep = "\t")  
}
\arguments{
  \item{base}{\code{base} a file name for a  matrix with two columns. 
    The first column is probe ids and the second one are the mappings to 
    SGD ids used by all the Yeast Genome data files. If \code{base} = "", 
    the whole genome will be mapped based on a data file that contains 
    mappings between all the ORFs and SGD ids} 
  \item{GEOAccNum}{\code{GEOAccNum} a character string for the accession
    number given by GEO for a yeast data set}
  \item{GEOUrl}{\code{GEOUrl} a character string for the url that
    contains a common CGI for all the GEO data. Currently it is
    \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?}}
  \item{geoCols}{\code{geoCols} a vector of integers for the coloumn
    numbers of the source file from GEO that maps yeast probe ids to ORF
    ids} 
  \item{yGenoUrl}{\code{yGenoUrl} a character string for the url that is
    a directory in Yeast Genom web site that contains directories for
    yeast annotation data. Currently it is
    \url{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/}}
  \item{baseURL}{see yGenoUrl}
  \item{yGenoNames}{\code{yGenoNames} a vector of character strings for
    the names of yeast annotation data. Each of the strings can be
    appended to yGenoUrl to make a complete url for a data file}
  \item{fileName}{a character string for the extension part of the
    source data file that can be used to target genes to SGD ids}
  \item{toKeep}{\code{toKeep} a list of vector of integers with numbers
    corresponding to column numbers of yeast genom data files that will
    kept when data files are processed. The length of toKeep must be the
    same as yGenoName (a vector for each file)}
  \item{colNames}{\code{colNames} a list of vectors of character strings
    for the names to be given to the columns to keep when processing the
    data. Again, the length of colNames must be the same as yGenoNames}
  \item{seps}{\code{seps} a vector of characters for the separators used
    by the data files included in yGenoNames}
  \item{sep}{singular version of seps}
  \item{by}{\code{by} a character string for the column that is common
    in all data files to be processed. The column will be used to merge
    separate data files}
  \item{probe2ORF}{\code{probe2ORF} a matrix with mappings of yease
    target genes to ORF ids that in turn can be mapped to SGD ids}
  \item{gos}{\code{gos} a vector of character strings for GO ids
    retrieved from Yeast Genome Project}
  \item{evis}{\code{evis} a vector of character string for the evidence
    code associated with go ids}
  \item{chr}{\code{chr} a vector of character strings for chromosome
    numbers} 
  \item{chrloc}{\code{chrloc} a vector of integers for chromosomal
    locations}
  \item{chrori}{\code{chrori} a vector of characters that can either be
    w or c that are used for strand of yeast chromosomes}
  \item{srcUrl}{\code{srcUrl} a character string for the url where
    source yeast genome data are stroed}
  \item{yGenoName}{\code{yGenoName} a character string for the yeast
    genome file name to be processed}
}
\details{
  To merge files, the system has to map the target genes in the base
  file to SGD ids and then use SGD ids to map traget genes to annotation
  data from different sources.

  \code{\link{formatGO}} adds leading 0s to goids when needed and then
  append the evidence code to the end of a goid following a "@".

  \code{\link{formatChrLoc}} assigns a + or - sing to \code{chrloc}
  depending on whether the corresponding \code{chrori} is w or c and
  then append \code{chr} to the end of \code{chrloc} following a "@".

  \code{\link{getGEOYeast}} gets yeast data from GEO for the columns
  specified. 
}
\value{
  \code{\link{yeastAnn}} returns a matrix with traget genes annotated by
  data from selected data columns in different data sources.

  \code{\link{getProbe2SGD}} returns a matrix with mappings between
  target genes and SGD ids.

  \code{\link{procYeastGeno}} returns a data matrix.

  \code{\link{formatGO}} returns a vector of character strings.

  \code{\link{formatChrLoc}} returns a vector of character strings.

  \code{\link{getGEOYeast}} returns a matrix with the number of columns
  specified. 
}
\references{\url{ftp://genome-ftp.stanford.edu}}
\author{Jianhua Zhang}

\seealso{}
\examples{
\dontrun{
yeastData <- yeastAnn(GEOAccNum = "GPL90")
}
}
\keyword{manip}
