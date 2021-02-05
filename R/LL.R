# Constructs an object of LL
#

LL <- function(srcUrl = getSrcUrl("LL"),
               parser = file.path(.path.package("AnnBuilder"), "scripts",
               "gbLLParser"), baseFile = "", built, fromWeb = TRUE){

    new("LL", srcUrl = srcUrl, parser = parser, baseFile = baseFile,
        built = ifelse(missing(built), getSrcBuilt("LL"), built),
        fromWeb = fromWeb)
}

# This function parses a source file based on a given set of
# parsing instrucitons passed as "parser" and "baseFile".
#
# Copyright 2002. J. Zhang all rights reserved.
#

fileMuncher <- function(outName, baseFile, dataFile, parser, isDir = FALSE){
    OS <- .Platform$OS.type
    perlName <- paste(tempfile("tempPerl"), "pl", sep=".")

    writePerl <- function(toWrite){
        write(toWrite, file = perlName, append = TRUE)
    }

    ## Convert \ => / for Windows
    if (OS == "windows") {
        outName <- chartr("\\", "/", outName)
        baseFile <- chartr("\\", "/", baseFile)
        dataFile <- chartr("\\", "/", dataFile)
    }

    if(!file.create(perlName))
        stop(paste("You do not have write permission to the ",
                   "directory for the Perl script createded", sep = ""))
    if(OS == "unix"){
        perlBin <- system("which perl", intern = TRUE)
        if(length(perlBin) > 1){
            stop("Perl is not available!")
        }
        writePerl(paste("#!", perlBin, "\n\n", sep = ""))
    }else if(OS == "windows"){
        writePerl("#!/usr/bin/perl -w\n\n")
    }

    if(baseFile != "" && !is.null(baseFile) && !is.na(baseFile)){
        statement <- paste("open(BASE, \"<", baseFile, "\") || die ",
                       "\"Can not open ", baseFile, "\";\n", sep = "")
        writePerl(statement)
    }

    if(dataFile != "" && !is.null(dataFile) && !is.na(dataFile)){
        if(isDir){
            statement <- paste("opendir(DIR, \"", dataFile, "\") || die ",
                           "\"Can not open director ",dataFile, "\";\n",
                           sep = "")
            writePerl(statement)
        }else{
            statement <- paste("open(DATA, \"<", dataFile, "\") || die ",
                       "\"Can not open ", dataFile, "\";\n", sep = "")
            writePerl(statement)
        }
    }

    statement <- paste("open(OUT, \">", outName, "\") || die ",
                       "\"Can not open ", outName, "\";\n\n", sep = "")
    writePerl(statement)
    if(isDir)
        writePerl(paste("$PATH = \"", dataFile, "\";", sep = ""))

    if(!is.null(parser))
        writePerl(readLines(parser))

    .callPerl(perlName, OS)
    unlink(perlName)
    return (outName)
}

.callPerl <- function(script, os){
    if(os == "unix"){
        system(paste("chmod +x", script))
        system(script)
    }else if(os == "windows"){
        script <- gsub("/", "\\\\", script)
        system(paste("perl", script))
    }else{
        stop(paste("Do not know who to run perl under ", os))
    }
}








