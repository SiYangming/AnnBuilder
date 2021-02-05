# This function creates descriptive statistics for a given data
# package.
#
# Copyright 2002, J. Zhang. All rights reserved.
#

getDPStats <- function(baseF, pkgName, pkgPath, saveList = TRUE,
                       isFile = TRUE){

    dateBuilt <- getDate(pkgName, pkgPath)
    keyMapped <- getProbeNum(pkgName, pkgPath, TRUE)

    if(is.null(baseF) || is.na(baseF) || baseF == ""){
        probeMapped <- NULL
        numMissMatch <- NULL
        srcNumProbes <- NULL
        numProbes <- NULL
        keyMissMatch <- NULL
        otherMapped <- keyMapped[["mapped"]]
    }else{
        probeBased <- paste(pkgName, getPBased(), sep = "")
        mappedKeys <- keyMapped[["mapped"]]
        otherMapped <- mappedKeys[!is.element(names(mappedKeys), probeBased)]
        if(isFile){
            srcNumProbes <- length(readLines(baseF))
        }else{
            srcNumProbes <- nrow(baseF)
        }
        numProbes <- keyMapped[["total"]]
        keyMissMatch <- matchProbes(baseF, pkgName, pkgPath,
                                    probeBased, isFile = isFile)
        numMissMatch <- NULL
        for(i in names(numProbes)){
            if(any(i == probeBased) && numProbes[i] != srcNumProbes){
                numMissMatch <- c(numMissMatch, i)
            }
        }
        if(is.null(numMissMatch)){
            numMissMatch <- "None"
        }

        probeMapped <- NULL
        for(i in names(mappedKeys)){
            if(any(i == probeBased)){
                tempNum <- as.numeric(mappedKeys[i])
                names(tempNum) <- i
                probeMapped <- c(probeMapped, tempNum)
            }
        }
    }
    statsList <- list(name = pkgName,
                      built = dateBuilt,
                      probeNum = srcNumProbes,
                      numMissMatch = numMissMatch,
                      probeMissMatch = keyMissMatch,
                      probeMapped = probeMapped,
                      otherMapped = otherMapped,
                      sessionInfo=sessionInfo())
    QCList <- formatABQCList(statsList)
    if(saveList){
        dataPath <- file.path(pkgPath, pkgName, "data")
        
        qcVar <- paste(pkgName, "QC", sep = "")
        qcFile <- file.path(dataPath, paste(qcVar, "rda", sep="."))
        assign(qcVar, QCList)
        save(list=qcVar, file=qcFile)

        qcDatVar <- paste(pkgName, "QCDATA", sep="")
        assign(qcDatVar, statsList)
        qcDatFile <- file.path(dataPath, paste(qcDatVar, "rda", sep="."))
        save(list=qcDatVar, file=qcDatFile)
        
        mapCountVar <- paste(pkgName, "MAPCOUNTS", sep = "")
        mapCountFile <- file.path(dataPath, paste(mapCountVar, "rda", sep="."))
        assign(mapCountVar, c(probeMapped, otherMapped))
        save(list=mapCountVar, file=mapCountFile)
    }
    return(statsList)
}

formatABQCList <- function(x){
    QCList <- paste("\n\nQuality control information for ", x$name,
                    "\nDate built:", x$built, "\n",
                    ifelse(is.null(x$probeNum), "",
                           paste("\nNumber of probes:", x$probeNum)),
                    ifelse(is.null(x$numMissMatch), "",
                           paste("\nProbe number missmatch:",
                                 paste(x$numMissMatch, sep = "",
                                       collapse = "; "))),
                    ifelse(is.null(x$probeMissMatch), "",
                           paste("\nProbe missmatch:",
                                 paste(x$probeMissMatch, sep = "",
                                       collapse = "; "))),
                    ifelse(is.null(x$probeMapped), "",
                           paste("\nMappings found for ",
                                 "probe based rda files: \n",
                                 paste("\t", names(x$probeMapped),
                                       "found", x$probeMapped,
                                       "of", x$probeNum, sep = " ",
                                       collapse = "\n"),
                                 sep = "", collapse ="\n")),
                    "\nMappings found for non-probe based rda files:\n",
                    paste("\t", names(x$otherMapped), "found",
                          x$otherMapped, sep = " ", collapse = "\n"),"\n\n")
    return(QCList)
}

getDate <- function(pkgName, pkgPath, fromDesc = TRUE){
    if(fromDesc){
        fileName <- file.path(pkgPath, pkgName, "DESCRIPTION")
        if(!file.exists(fileName))
            stop(paste("DESCRIPTION file was not found in",
                       file.path(pkgPath, pkgName), "!"))

        tempFile <- readLines(fileName)

        for (i in tempFile){
            if(regexpr("Created:", i) > 0){
                return(gsub("^Created:\\(*\\)", "\\", i))
            }
        }
        return("Date built not found")
    }else{
        return(date())
    }
}

getProbeNum <- function(pkgName, pkgPath, noNA = FALSE){

    toReturn <- NULL
    dataFolder <- file.path(pkgPath, pkgName, "data")
    rdaFiles <- getDirContent(dataFolder)
    nums <- sapply(rdaFiles, countMapping, noNA = FALSE)
    mapped <- sapply(rdaFiles, countMapping, noNA = TRUE)
    names(nums) <- gsub("\\..*", "", basename(names(nums)))
    names(mapped) <- gsub("\\..*", "", basename(names(nums)))

    return(list(total = nums, mapped = mapped))
}

countMapping <- function(rdaName, noNA = FALSE){
    #print(paste("Processing", rdaName))
    tempEnv <- new.env(hash = TRUE, parent = NULL)
    rda <- gsub("\\.rda", "", basename(rdaName))
    load(rdaName, tempEnv)

    tempData <- tempEnv[[rda]]
    if(is.environment(tempData)){
        tempData <- unlist(eapply(tempEnv[[rda]],
                                  FUN = function(x) !is.na(x)[1]))
    #    tempData <- as.list(tempData)
    }else{
        tempData <- unlist(sapply(tempData, FUN = function(x) !is.na(x)))
    }
    if(!noNA){
        return(length(tempData))
    }else{
        #fun <- function(x) !is.na(x)[1]
        #temp <- unlist(sapply(tempData, fun))
        return(length(tempData[tempData]))
    }
}

matchProbes <- function(baseF, pkgName, pkgPath, toMatch, isFile = TRUE){
    tempEnv <- new.env(hash = TRUE, parent = NULL)
    toReturn <- NULL

    if(isFile){
        srcProbes <- read.delim(baseF, header = FALSE, sep = "\t")
    }else{
        srcProbes <- baseF
    }
    dataFolder <- file.path(pkgPath, pkgName, "data")
    rdaFiles <- list.files(dataFolder)

    for(i in rdaFiles){
        if(regexpr(".rda", i) > 0){
            name <- gsub(paste(pkgName,"(.*).rda", sep = ""),
                         "\\1", i)
            load(file.path(dataFolder, i), tempEnv)
            if(any(toMatch == name)){
               if(any(!is.element(srcProbes[,1],
                         names(as.list(tempEnv[[gsub(".rda", "", i)]]))))){
                    toReturn <- c(toReturn, i)
                }
            }
        }
    }
    return(ifelse(is.null(toReturn), "None", toReturn))
}

getPBased <- function(){

    info <- unlist(as.list(AnnInfo))
    pBased <- info[grep("pbased", names(info))]

    return(gsub(".pbased", "", names(pBased[pBased == "Y"])))
}
