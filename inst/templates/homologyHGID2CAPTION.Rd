\name{homologyHGID2CAPTION}
\alias{homologyHGID2CAPTION}
\title{Mappings between NCBI's Internal HomoloGeneIDs (HGID) to
  descriptive captions of the HGIDs}   
\description{
  homologyHGID2CAPTION provides descriptive captions of HomoloGeneIDs
  (HGID) used by NCBI to represent sequences in various organisms. A
  HGID represents a sequence that is conserved in a group of organisms.   
}
\details{
  Each HGID has one descriptive caption but a given descriptive caption
  may include multiple HGIDs. This is because homologous relations may be
  found in different groups of organisms or different genes in the same
  group of organisms. For example, HGID 37237 and 37238 are both
  described as 'Gene conserved in Mammalia'. However, the former represents
  a sequence found in H.sapiens (LOC402318) and R.norvegicus (LOC362421)
  while the latter in H.sapiens (LOC402318), H.sapien2 (LOC391370), and
  R.norvegicus (LOC362421) (based on data downloaded from HomoloGene on
  08/11/2004).  
 
  The Entrez Gene identifiers (across organisms) for the sequence represented by a
  given HGID can be obtained using the environment homologyHGID2LL.  
  
  Mappings contained were based on data provided by HomoloGene.
  
  #HGSOURCE#

  Packagebuilt #DATE#
}
\references{
\url{http://www.ncbi.nih.gov/entrez/query.fcgi?db=homologene}
}
\examples{
    require("annotate") || stop("annotate unavailable")
    xx <- as.list(homologyHGID2CAPTION)
        if(length(xx) > 5){
                # Get the caption for one HGID
                xx[2]
                # Get the caption for a few HGID
                if(length(xx) >= 3){
                        xx[2:4]
                }
        }
}
\keyword{datasets} 

