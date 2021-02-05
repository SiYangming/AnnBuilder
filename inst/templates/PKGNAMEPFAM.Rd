\name{#PKGNAME#PFAM}
\alias{#PKGNAME#PFAM}
\title{Map Manufacturer IDs to Pfam IDs}
\description{
  #PKGNAME#PFAM is an R environment that provides mappings between a
 manufacturer identifier and the associated Pfam identifiers.
}
\details{
  Each manufacturer identifier maps to a named vector of Pfam identifiers.  The name
  for each Pfam identifier is the IPI accession numbe where this Pfam identifier is found.

  If the Pfam is a named NA, it means that the associated Entrez Gene id
  of this manufacturer identifier is found in an IPI entry of the IPI database,
  but there is no Pfam identifier in the entry.

  If the Pfam is a non-named NA, it means that the associated Entrez Gene id
  of this manufacturer identifier is not found in any IPI entry of the IPI database.
}
\references{
  Kersey P. J., Duarte J., Williams A., Karavidopoulou Y., Birney E., Apweiler R.
  The International Protein Index: An integrated database for proteomics experiments.
  Proteomics 4(7): 1985-1988 (2004)
  \url{http://www.ebi.ac.uk/IPI/IPIhelp.html}
}
\examples{
        # Convert to a list
        xxx <- as.list(#PKGNAME#PFAM)
        # randomly display 10 probes
        sample(xxx, 10)
}
\keyword{datasets}

