\name{pfam#FROM#2#TO#}
\alias{pfam#FROM#2#TO#}
\title{Mappings from #FROMNAME# to #TONAME#}
\description{
  Mapping function from #FROMNAME# (#FROM#) to the associated #TONAME# (#TO#)
}
\usage{
pfam#FROM#2#TO#(#FOR#=NULL)
}
\arguments{
  \item{#FOR#}{A character vector of #FROMNAME# or NULL (default)}
}
\value{
  A list of #TONAME#.
}
\details{
  If \code{#FOR#} is a vector of #FROMNAME#,  pfam#FROM#2#TO# will
  return a list of associated #TONAME# for the given \code{#FOR#}.  If
  \code{#FOR#} is NULL,  pfam#FROM#2#TO#(#FOR#=NULL) (or pfam#FROM#2#TO#())
  will return a list of associated #TONAME# for all the #FROMNAME# in
  the database.
  Extra information:  #EXTRA#
  More Details: #DETAILINFO#
}
\references{\url{http://www.sanger.ac.uk/Software/Pfam/} and
  \url{ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current\_release/userman.txt}}
\author{Ting-Yuan Liu}
\examples{
\dontrun{
    #FROM#2#TO# <- pfam#FROM#2#TO#()
    head(#FROM#2#TO#)
    pfam#FROM#2#TO#(#FOR#=sample(names(#FROM#2#TO#), 3))
}
}
\keyword{datasets}
