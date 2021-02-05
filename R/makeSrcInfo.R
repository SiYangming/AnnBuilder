# This function creates an R environment in .GlobalEnv to hold source
# data information that will be used later by other functions. The
# environment object contains key-value pairs with
# keys being data names such as LOCUSID, ACCNUM and values being a
# list containing a short (sort diescription), long (long
# description), and source (data source) element.
#
# Copyright 2002 J. Zhang, all right reserved.
#

makeSrcInfo <- function(srcFile = ""){
    path <- .path.package("AnnBuilder")

    if(is.null(srcFile) || is.na(srcFile) || srcFile == ""){
        srcFile <- file.path(path, "data", "AnnInfo")
    }

    info <- matrix(scan(srcFile, what = "", sep = "\t", quiet = TRUE),
                   ncol = 5, byrow = TRUE)
    temp<- new.env(hash = TRUE, parent = NULL)
    for(i in 1:nrow(info)){
        assign(info[i,][1], list(short = info[i,][2],
                                 long = info[i,][3],
                                 src = info[i,][4],
                                 pbased = info[i,][5]),
               env = temp)
    }
    assign("AnnInfo", temp, env = .GlobalEnv)
}

getAllSrc <- function(){
    return(c("LL", "UG", "GO", "KEGG", "GP", "YG", "GEO"))
}
