.annBuilderParseSourceURLs <- function() {
    dataFile <- "sourceURLs.txt"
    cols <- rep("character", 2)
    df <- read.table(dataFile, header=TRUE, colClasses=cols)
    sourceURLs <- as.list(df$url)
    names(sourceURLs) <- df$name
    sourceURLs
}
sourceURLs <- .annBuilderParseSourceURLs()
rm(.annBuilderParseSourceURLs)
