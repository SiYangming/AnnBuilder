# Building a data package containing location data for genes on
# chromosomes
chrLocPkgBuilder <- function(pkgName = "humanCHRLOC", pkgPath, version,
                             author, organism = "Homo sapiens")
    {
    url = getSrcUrl("gp", organism)
    repList <- getRepList("CHRLOC",
                          list(gp = GP(srcUrl = url, organism = organism,
                               fromWeb = TRUE)))
    repList[["PKGNAME"]] <- pkgName

    chrLoc <- getChroLocation(url, raw = TRUE)
    suffix <- gsub(".*_random", "Unconfident", chrLoc[,2])
    suffix[suffix != "Unconfident" | is.na(suffix)] <- "Confident"
    chrLoc <- cbind(chrLoc[,1], gsub("_random", "", chrLoc[,2]),
                    chrLoc[,3], paste(chrLoc[, 4], suffix, sep = "@"),
                    paste(chrLoc[, 5], suffix, sep = "@"))
    # Put signs for - strands
    chrLoc[, 4] <- gsub("\\+", "", paste(chrLoc[,3], chrLoc[,4], sep = ""))
    chrLoc[, 5] <- gsub("\\+", "", paste(chrLoc[,3], chrLoc[,5], sep = ""))
#    chrLength <- findChrLength(organism)
    makeSrcInfo()
    createEmptyDPkg(pkgName, pkgPath)
    # Write the mappings between LL and chromosome number
    ll2Chr <- mergeRowByKey(chrLoc[, 1:2])
    saveMat(ll2Chr[, 1:2], fun = splitEntry, pkgName = pkgName,
                pkgPath = pkgPath, envName = "LOCUSID2CHR")

    for(i in getChrNum(unique(chrLoc[,2]))){
        chr4Org <- chrLoc[chrLoc[,2] == i,]
        # Remove <NA>s introduced by the above manipulation
        chr4Org <- chr4Org[!is.na(chr4Org[,1]), ]
        # For -strand start should be end and end should be start
        # chr4Org <- rbind(chr4Org[chr4Org[, 3] == "+", c(1, 4, 5)],
        #                 chr4Org[chr4Org[, 3] == "-", c(1, 5, 4)] )
        chr4Org <- chr4Org[, c(1,4,5)]
        # Collapse values for LocusLink ids
        chr4Org <- mergeRowByKey(chr4Org)
        fun <- function(x) twoStepSplit(x, asNumeric = TRUE)
        saveMat(chr4Org[, 1:2], fun = fun, pkgName = pkgName,
                pkgPath = pkgPath,
                envName = paste(i, "START", sep = ""))
        #env <- new.env(hash = TRUE, parent = NULL)
        #multiassign(chr4Org[, 1], lapply(chr4Org[, 2],
        #                               twoStepSplit, asNumeric = TRUE), env)
        #assign(paste(pkgName, i, "START", sep = ""), env)
        #save(list = paste(pkgName, i, "START", sep = ""),
        #     file =  file.path(pkgPath, pkgName, "data",
        #     paste(pkgName, i, "START.rda", sep = "")))
        saveMat(chr4Org[, 1:3], fun = fun, pkgName = pkgName,
                pkgPath = pkgPath,
                envName = paste(i, "END", sep = ""))
    }
    saveCytoband(pkgName, pkgPath, organism, url, ext = "cytoBand.txt.gz")

    for(i in getChroms4Org(organism)){
        repList[["CHROM"]] <- as.character(i)
        copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "PKGNAMEEND.Rd"),
                       file.path(pkgPath, pkgName, "man",
                                 paste(pkgName, i, "END.Rd",
                                                sep = "")), repList,   "#")
        copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "PKGNAMESTART.Rd"),
                       file.path(pkgPath, pkgName, "man",
                                 paste(pkgName, i, "START.Rd",
                                                sep = "")), repList,   "#")
    }
    copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "PKGNAMECYTOLOC.Rd"),
                       file.path(pkgPath, pkgName, "man",
                                 paste(pkgName,"CYTOLOC.Rd",
                                                sep = "")), repList,   "#")
    copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "PKGNAMELOCUSID2CHR.Rd"),
                       file.path(pkgPath, pkgName, "man",
                                 paste(pkgName,"LOCUSID2CHR.Rd",
                                                sep = "")), repList, "#")
    writeDescription(pkgName, pkgPath, version, author)
    #copyTemplates(repList, pattern, pkgName, pkgPath)
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates",
                             "GENERAL.Rd"), file.path(pkgPath, pkgName,
                             "man", paste(pkgName, ".Rd", sep = "")),
                              list(PKGNAME = pkgName, SOURCENBUILT =
                                   paste("\n\t", repList, sep = "",
                                         collapse = "\n")), "#")
    writeZZZ(pkgPath, pkgName)
    writeFun(pkgPath, pkgName)
    getDPStats("", pkgName, pkgPath)
    copySubstitute(file.path(.path.package("AnnBuilder"),
                        "templates", "PKGNAMEQC.Rd"),
                        file.path(pkgPath, pkgName, "man",
                        paste(pkgName, "QC.Rd", sep = "")),
                        list(PKGNAME = pkgName), "#")
    makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
}
# Chromosome number may have _random extensions for genes not too sure
# about their locations. Remove them
getChrNum <- function(chr){
    options(warn = -1)
    keep <- chr[!is.na(as.numeric(chr))]
    options(warn = 1)
    keep <- c(keep, chr[is.element(chr, c("X", "Y"))])
    return(keep)
}

saveCytoband <- function(pkgName, pkgPath, organism, url,
                         ext = "cytoBand.txt.gz"){
    cytoband <- getCytoLoc(organism, srcUrl = paste(url, ext, sep = ""))
    saveList(cytoband, pkgName, pkgPath, "CYTOLOC")
    writeManPage(pkgName, pkgPath, "CYTOLOC", organism = organism,
                 src = "gp", isEnv = FALSE)
}

getChroms4Org <- function(organism){
    switch(toupper(organism),
           "HOMO SAPIENS" = return(c(as.character(1:22), "X", "Y")),
           "MUS MUSCULUS" = return(c(as.character(1:19), "X", "Y")),
           "RATTUS NORVEGICUS" = return(c(as.character(1:20), "X", "Y")),
           stop("Unknown organism name"))
}

