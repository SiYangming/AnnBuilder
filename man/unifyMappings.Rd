\name{unifyMappings}
\alias{unifyMappings}

\title{A function to unify mapping result from different sources}
\description{
  Given a base file and mappings from different sources, this function
  resolves the differences among sources in mapping results using a
  voting scheme and derives unified mapping results for targets in the
  base file
}
\usage{
unifyMappings(base, eg, ug, otherSrc)
}

\arguments{
  \item{base}{\code{base} a matrix with two columns. The first column
    contains the target items (genes) to be mapped and the second the
    know mappings of the target to GenBank accession numbers or UniGene ids}
  \item{eg}{\code{eg} an object of class EG}
  \item{ug}{\code{ug} an object of class UG}
  \item{otherSrc}{\code{otherSrc} a vector of character strings for
    names of files that also contain mappings of the target genes in
    base. The files are assumed to have two columns with the first one
    being target genes and second one being the desired mappings}
}
\details{
  eg and ug have methods to parse the data from LocusLink and UniGene to
  obtain mappings to target genes in base. Correct source urls
  and parsers are needed to obtain the desired mappings 
}
\value{
  The function returns a matrix with four columns. The first two are the
  same as the columns of base, the third are unified mappings, and forth
  are statistics of the agreement among sources.
}

\author{Jianhua Zhang}
\note{This function is part of the Bioconductor project at Dana-Farber
  Cancer Institute to provide Bioinformatics functionalities through R}

\seealso{\code{\link{EG}}, \code{\link{UG}}}
\examples{
\dontrun{
myDir <- file.path(.path.package("AnnBuilder"), "temp")
geneNMap <- matrix(c("32468_f_at", "D90278", "32469_at", "L00693",
                   "32481_at", "AL031663", "33825_at", " X68733",
                   "35730_at", "X03350", "36512_at", "L32179",
                   "38912_at", "D90042", "38936_at", "M16652",
                   "39368_at", "AL031668"), ncol = 2, byrow = TRUE)
colnames(geneNMap) <- c("PROBE", "ACCNUM")
write.table(geneNMap, file = file.path(myDir, "geneNMap"), sep = "\t",
quote = FALSE, row.names = FALSE, col.names = FALSE)

temp <- matrix(c("32468_f_at", NA, "32469_at", "2",
                   "32481_at", NA, "33825_at", " 9",
                   "35730_at", "1576", "36512_at", NA,
                   "38912_at", "10", "38936_at", NA,
                   "39368_at", NA), ncol = 2, byrow = TRUE)
temp
write.table(temp, file = file.path(myDir, "srcone"), sep = "\t",
quote = FALSE, row.names = FALSE, col.names = FALSE)
temp <- matrix(c("32468_f_at", NA, "32469_at", NA,
                   "32481_at", "7051", "33825_at", NA,
                   "35730_at", NA, "36512_at", "1084",
                   "38912_at", NA, "38936_at", NA,
                   "39368_at", "89"), ncol = 2, byrow = TRUE)
temp
write.table(temp, file = file.path(myDir, "srctwo"), sep = "\t",
quote = FALSE, row.names = FALSE, col.names = FALSE)
otherMapping <- c(srcone = file.path(myDir, "srcone"),
srctwo = file.path(myDir, "srctwo"))

baseFile <-  file.path(myDir, "geneNMap")
egParser <- file.path(.path.package("AnnBuilder"), "scripts", "gbLLParser")
ugParser <- file.path(.path.package("AnnBuilder"), "scripts", "gbUGParser")
#if(.Platform$OS.type == "unix"){
    egUrl <-  "http://www.bioconductor.org/datafiles/wwwsources"
    ugUrl <-  "http://www.bioconductor.org/datafiles/wwwsources/Ths.data.gz"
    fromWeb = TRUE
#}else{
#    egUrl <- file.path(.path.package("AnnBuilder"), "data", "Tll_tmpl")
#    ugUrl <- file.path(.path.package("AnnBuilder"), "data", "Ths.data")
#    fromWeb = FALSE
#}
eg <- EG(srcUrl = egUrl, parser = egParser, baseFile = baseFile,
accession = "Tll_tmpl.gz")
ug <- UG(srcUrl = ugUrl, parser = ugParser, baseFile = baseFile,
organism = "Homo sapiens") 
# Only works interactively
    unified <- unifyMappings(base =  geneNMap, eg = eg, ug = ug,
               otherSrc = otherMapping)
    read.table(unified, sep = "\t", header = FALSE)

    unlink(c(file.path(myDir, "geneNMap"), file.path(myDir, "srcone"),
    file.path(myDir, "srctwo"), unified))
}
}
\keyword{manip}



