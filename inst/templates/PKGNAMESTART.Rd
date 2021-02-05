\name{#PKGNAME##CHROM#START}
\alias{#PKGNAME##CHROM#START}
\title{An annotation data file for transciption starting locations of genes
on chromosome #CHROM#}
\description{
  #PKGNAME##CHROM#START maps LocusLink ids to the transciption starting
  locations of genes on chromosome number #CHROM# corresponding to the
  LocusLink ids 
}
\details{
   This is an environment object containing key and value pairs. Keys
   are LocusLink ids and values are the transciption starting location for
   genes. The starting positions for genes on both the sense and antisense
   strand are number of base pairs measured from the p (5' end of the
   sense strand) to q (3' end of the sense strand) arms. Values for the
   antisense strand have a leading "-" sign (e. g. -1234567).

   Values for keys are named vectors of length 1 derived from the
   smallest starting value found for the LocusLink id. Names of chromosome
   location values can be "Confident" when the gene can be confidently
   placed on a chromosome and "Unconfident" otherwise (denoted by
   "\_random" in UCSC data).

   Mappings were derived from the following public data sources:

   #GPSOURCE#

   Package built: #DATE#
}
\references{
\url{http://www.genome.ucsc.edu/goldenPath/hg16/database/}
}
\examples{
	require("annotate") || stop("annotate unavailable")
	xx <- as.list(#PKGNAME##CHROM#START)
	if(length(xx) > 0){
		# Get the value of the first key
		xx[1]
		# Get the values for a few keys
		if(length(xx) >= 3){
		    xx[1:3]
		}
	}
}
\keyword{datasets}

