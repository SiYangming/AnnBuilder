# Constructors an object of GEO
#

GEO <- function(srcUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
                built){

    new("GEO", srcUrl = srcUrl, built = ifelse(missing(built),
                                getSrcBuilt("GEO"), built))
}


# Query the GEO database. srcUrl gives the url for a common CGI scrips
# and GEOAccNum is the GEO accession number representing a file in the
# database
queryGEO <- function(GEOObj, GEOAccNum){

    srcURL <- paste(srcUrl(GEOObj), "acc=", GEOAccNum,
                     "&view=data&form=text&targ=self", sep = "")

    conn <- url(srcURL, open = "r")
    temp <- readLines(conn)
    close(conn)
    # Remove the header lines that come with the file
    temp <- temp[grep("\t", temp)]
    # Add NAs to lines with no value for the last column
    temp <- strsplit(gsub("\t$", "\tNA", temp), "\t")
    # Convert to a matrix
    temp <- t(sapply(temp, unlist))
    # The first row is for column name. Remove it.
    colnames(temp) <- temp[1,]
    return(temp[-1,])
}
