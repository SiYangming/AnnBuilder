# Constructors an object of GEO
#

YG <- function(srcUrl =
                  "ftp://genome-ftp.stanford.edu/pub/yeast/data_download/",
                      parser = "", baseFile = "", built, fromWeb = TRUE){

    new("YG", srcUrl = srcUrl, parser = parser, baseFile = baseFile,
        built = ifelse(missing(built), getSrcBuilt("YG"), built),
        fromWeb = fromWeb)
}

# Extracts the speficed columns from yeast genomic data. srcUrl is the base
# url for the ftp site, extenName is the extension for a target file
# name, cols2Keep is a vector for the number of data columns to keep,
# and sep is the delimiter used by the file.

getYeastData <- function(url, extenName, cols2Keep, sep){
    extenFilename <- loadFromUrl(paste(url, extenName, sep=""))
    conn <- file(extenFilename, open="r")
    ##conn <- url(paste(url, extenName, sep = ""), open = "r")
    skip <- 0
    if(length(grep(".*gene_association.sgd", extenName)) > 0){
        comment <- "!"
    }else{
        comment <- ""
    }
    options(show.error.messages = FALSE)
    # Try read table first
    data <- try(read.table(conn, header = FALSE, sep = sep, quote = "",
                           skip = skip, comment.char = comment,
                           strip.white  = TRUE))
    options(show.error.messages = TRUE)
    close(conn)
    # Has to use readLines as the source data are not well structured
    # and columns may miss if data = "try-error"
    if(inherits(data, "try-error")){# | length(cols2Keep) > dim(data)[2]){
        data <- readBadData(paste(url, extenName, sep = ""), sep)
    }
    data <- as.matrix(data)
    return(data[, cols2Keep])
}


# Some of the yest data files are badly structured (e.g. missing
# columns) and have to be processed this way
readBadData <- function(url, sep){
     conn <- url(url, open = "r")
     data <- readLines(conn)
     close(conn)
     if(length(data) <= 10){
         colNum <- findNumCol(data, sep = sep)
     }else{
         colNum <- findNumCol(data[1:10], sep = sep)
     }

     data <- as.matrix(sapply(data, strsplit, sep))
     lengths <- as.vector(sapply(data, length))
     data <- unlist(data[lengths == colNum], use.names = FALSE)

     return(matrix(data, ncol = colNum, byrow = TRUE))
}

# Figures out the number of columns based on a section of a file
findNumCol <- function(fewLines, sep){
    lengths <- sapply(sapply(fewLines, strsplit, sep), length)
    return(max(lengths[duplicated(lengths)]))
}

