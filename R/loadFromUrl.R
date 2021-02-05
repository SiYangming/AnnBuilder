# This function tries to download a file from a given target source
# and unzip the file if needed.
#
# Copyright 2002, J. Zhang. All rights reserved.
#
loadFromUrl <- function(srcUrl, destDir = "", verbose=FALSE){
    if( verbose )
        cat("loading from URL:", srcUrl, "\n")

    if(destDir == ""){
        destDir <- tempdir()
    }

    fileName <- file.path(destDir, paste(basename(tempfile()),
                             gsub("\\.*.*/(\\.*.*)", "\\1", srcUrl), sep = ""))
    # Make sure writing is permited
    options(show.error.messages = FALSE)
    tryMe <- try(file(fileName, "w"))
    options(show.error.messages = TRUE)
    if(inherits(tryMe, "try-error")){
        stop(paste("You do not have write permission to",
                   "the direcory specified!"))
    }
    close(tryMe)
    options(show.error.messages = FALSE)
    if(.Platform$OS.type=="unix"){
      tryMe <- try(download.file(srcUrl, fileName,
                      method = "internal", quiet = TRUE))
    }
    if(.Platform$OS.type=="windows"){
      tryMe <- try(download.file(srcUrl, fileName,
                                 method="internal", mode="wb",
                                 quiet=TRUE))
    }
    options(show.error.messages = TRUE)
    if(inherits(tryMe, "try-error")){
        stop(paste("URL", srcUrl, "is incorrect or the target",
                   "site is not responding!"))
    }

    # Unzip if necessary
    if(regexpr("\\.*.gz", fileName) > 0){
        if(file.exists(gsub("(\\.*.*).gz", "\\1", fileName)))
            unlink(gsub("(\\.*.*).gz", "\\1", fileName))
        unzipFile(fileName, destDir, TRUE)
    }else{
        if(regexpr("\\.*.zip", fileName) > 0){
            unzipFile(fileName, destDir, FALSE)
        }
    }
    return (gsub("(.gz)|(.zip)", "\\", fileName))
}

validateUrl <- function(srcUrl){
    options(show.error.messages = FALSE)
    con <- try(url(srcUrl, "r"))
    options(show.error.messages = TRUE)
    if(inherits(con, "try-error")){
        stop(paste("URL", srcUrl, "is incorrect or the target",
                   "site is not responding!"))
    }else{
        close(con)
    }
}

unzipFile <- function(fileName, where =
                      file.path(.path.package("AnnBuilder"), "data"),
                      isgz = FALSE){
    curDir <- getwd()
    setwd(where)
    if(.Platform$OS.type == "unix"){
        if(isgz){
            system(paste("gunzip -q", fileName), intern=TRUE)
        }else{
            if(!is.null(getOption("unzip"))){
                system(paste(getOption("unzip"), "-q", basename(fileName)), intern=TRUE)
            }else{
                stop("Can not find unzip in the path")
            }
        }
    }else if(.Platform$OS.type == "windows"){
        if(isgz){
            tst <- try(system(paste("gzip -d",fileName), intern=TRUE))
            options(show.error.messages = TRUE)
            if(inherits(tst, "try-error")){
                stop(paste("Gzip is either not installed or not in your path.\n",
                           "You may need to install the Rtools available at:\n",
                           " http://www.stats.ox.ac.uk/pub/Rtools/"))
              }
        }else{
            zip.unpack(fileName, getwd())
        }
        unlink(fileName)
    }
    setwd(curDir)
}

