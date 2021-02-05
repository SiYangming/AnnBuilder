GOPkgBuilder <- function(pkgName, pkgPath, filename, version, author,
                         lazyLoad=TRUE)
{
    sourceURLs <- getOption("AnnBuilderSourceUrls")
    eg = EG(srcUrl=sourceURLs[["EG"]])
    go = GO(srcUrl=sourceURLs[["GO"]])
    repList = getRepList("GO", list(eg=eg, go=go))
    createEmptyDPkg(pkgName, pkgPath, force=TRUE)
    
    cat("creating main environments...")
    godata = GOXMLParser(filename)
    print(paste(class(godata[["TERM"]]), length(godata[["TERM"]]), sep=" "))
    for (dataName in names(godata)) {
        cat("saving", dataName, "\n")
        saveEnv(godata[[dataName]], pkgName=pkgName, pkgPath=pkgPath,
                envName=dataName)
    }
    cat("DONE\n")
    
    egDat = getEG(godata[["OBSOLETE"]])
    gl = getGOLOCUSID(egDat)
    saveEnv(gl, pkgName = pkgName, pkgPath = pkgPath, envName = "LOCUSID")

    gl2g = createGoLocusId2Go(egDat, godata[["TERM"]])
    saveEnv(gl2g, pkgName = pkgName, pkgPath = pkgPath, envName = "LOCUSID2GO")

    offSpringEnv = new.env(parent=NULL)
    l2e(as.list(godata[["BPOFFSPRING"]]), offSpringEnv)
    l2e(as.list(godata[["MFOFFSPRING"]]), offSpringEnv)
    l2e(as.list(godata[["CCOFFSPRING"]]), offSpringEnv)
    
    saveEnv(getGOALLLOCUSID(gl, offSpringEnv), pkgName=pkgName,
            pkgPath=pkgPath, envName="ALLLOCUSID")
    
    writeDescription(pkgName, pkgPath, version, author)
    writeZZZ(pkgPath, pkgName)
    ## Write the quality control data
    repList[["PKGNAME"]] <- "GO"
    writeDocs("", pkgName, pkgPath, version, author, repList, "GO")
    if (lazyLoad) {
        makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
    }
}


getEG <- function(goObsolete)
{
    NUM_COLS <- 7
    keepCols <- c(2, 3, 4)
    
    sourceURLs <- getOption("AnnBuilderSourceUrls")
    src <- loadFromUrl(file.path(sourceURLs[["EG"]], "gene2go.gz"))
    input = read.table(src, na.strings="-", sep="\t", quote="",
      comment.char="#", colClass=rep("character", NUM_COLS), header=FALSE)
    ## include 7 species human, mouse, rat, drosophila, zebrafish, Saccharomyces
    ## cerevisiae and Caenorhabditis elegans
    wantedSpecies <- c("9606", "10090", "10116", "6239", "4932",
                       "7955", "7227")
    dat = input[(input[,1] %in% wantedSpecies), keepCols]
    names(dat) <- c("egId", "goId", "evi")

    ## Remove rows containing obsolete GOIds
    obsoleteFound <- dat$goId %in% ls(goObsolete)
    dat <- dat[!obsoleteFound, ]
    if (sum(obsoleteFound) > 0)
      warning("Found ", sum(obsoleteFound), " obsolete GOIds in gene2go")
    return(dat)
}

getGOLOCUSID <- function(dat)
{
    vec = dat[,1]
    names(vec) = dat[,3]
    mylist = split(vec,dat[,2])
    myEnv = new.env(parent=NULL)
    l2e(mylist, myEnv)
    return(myEnv)
}


createGoLocusId2Go <- function(gene2goDf, goTermEnv)
{
    ## Return a three element list with names BP, CC, MF.  Each element is an
    ## environment mapping Entrez Gene Id to a vector of associated GO Ids named
    ## with evidence code.
    expectedNames <- c("egId", "goId", "evi")
    if ((ncol(gene2goDf) != length(expectedNames)) ||
        (names(gene2goDf) != expectedNames))
      stop("Expecting ", length(expectedNames), " columns with names ",
           expectedNames, ".  Got: ", names(gene2goDf))
    
    ## filter out missing GOIds
    knownGoIds <- gene2goDf$goId %in% ls(goTermEnv)
    if (sum(!knownGoIds) > 0) {
        badIds <- paste(gene2goDf$goId[!knownGoIds], collapse=", ")
        warning("Dropping records from gene2go with GO Ids not found in GOTERM environment:\n",
                badIds)
    }
    gene2goDf <- gene2goDf[knownGoIds, ]
    goCat <- sapply(mget(gene2goDf$goId, goTermEnv), function(x) Ontology(x))
    splitData <- split(gene2goDf, goCat)
    rm(gene2goDf)

    nameAndSplitToEnv <- function(df) {
        goIds <- df$goId
        names(goIds) <- df$evi
        egIdSplit <- split(goIds, df$egId)
        env <- new.env(hash=TRUE, parent=NULL)
        l2e(egIdSplit, env)
    }
    splitData <- lapply(splitData, nameAndSplitToEnv)
    eg2goEnv <- new.env(hash=TRUE, parent=NULL)

    makeListStruct <- function(x, ont) {
        result <- list()
        for (i in 1:length(x)) 
          result[[i]] <- list(GOID=x[i], Evidence=names(x)[i], Ontology=ont)
        names(result) <- x
        result
    }

    for (cat in names(splitData)) {
        for (egId in ls(splitData[[cat]])) {
            ans <- makeListStruct(splitData[[cat]][[egId]], cat)
            eg2goEnv[[egId]] <- c(eg2goEnv[[egId]], ans)
        }
    }
    eg2goEnv
}


getGOALLLOCUSID <- function(dat, offspringEnv)
{
    ## need to merge offspring
    ## take goVec from offspring
    ## mget on GOLOCUSID
    
    myEnv = new.env(parent=NULL)
    goVec <- unlist(names(as.list(offspringEnv)))
    print(paste("length",length(goVec),sep=" ")) 
    for (i in 1:length(goVec))
      {
          mykey = as.character(goVec[i])
          ##print(mykey)
          myList = get(mykey, offspringEnv)
          myList = c(mykey,myList)
          if (!is.na(myList[[1]]))
            {
    		nextList = unlist(unname(mget(myList, dat, ifnotfound=NA)))
                ##print(nextList)
    		w = paste(nextList,names(nextList),sep="") 
    		if (length(nextList[!duplicated(w) & !is.na(nextList)]) == 0)
                  {
                      myEnv[[mykey]] <- NA
                  }
    		else
                  {
                      myEnv[[mykey]] <- nextList[!duplicated(w) & !is.na(nextList)]
                  }
            }
      }
    return (myEnv)
}


saveEnv <- function(envir, pkgName, pkgPath, envName) {
    env <- new.env(hash = TRUE, parent = NULL)
    copyEnv(envir, env)
    lockEnvironment(env, bindings = TRUE)
    assign(paste(pkgName, envName, sep = ""), env)
    save(list = paste(pkgName, envName, sep = ""), file = file.path(pkgPath,
                                                     pkgName, "data", paste(pkgName, envName, ".rda", sep = "")))
}


createEmptyDPkg <- function(pkgName, pkgPath,
                            folders = c("man", "R", "data"), force = TRUE) {
    if (file.exists(file.path(pkgPath, pkgName))) {
        if (!force) {
            stop(paste("Package", pkgName, "already exists"))
        }else{
            unlink(file.path(pkgPath, pkgName), TRUE)
        }
    }
    dir.create(file.path(pkgPath, pkgName))
    for (i in folders) {
        dir.create(file.path(pkgPath, pkgName, i))
    }
}


writeDocs <- function(baseName, pkgName, pkgPath, version, author,
                      repList, pattern, isFile = TRUE) {
    ## Write man pages and other files
    writeDescription(pkgName, pkgPath, version, author)
    copyTemplates(repList, pattern, pkgName, pkgPath)
    
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates",
                             "GENERAL.Rd"), file.path(pkgPath, pkgName,
                                                      "man", paste(pkgName, ".Rd", sep = "")),
                   list(PKGNAME = pkgName, SOURCENBUILT =
                        paste("\n\t", repList, sep = "",
                              collapse = "\n")), "#")
    writeZZZ(pkgPath, pkgName)
    writeFun(pkgPath, pkgName)
    getDPStats(baseName, pkgName, pkgPath, isFile)
    copySubstitute(file.path(.path.package("AnnBuilder"),
                             "templates", "PKGNAMEQC.Rd"),
                   file.path(pkgPath, pkgName, "man",
                             paste(pkgName, "QC.Rd", sep = "")),
                   list(PKGNAME = pkgName), "#")
}


copyTemplates <- function(repList, pattern, pkgName, pkgPath,
                          replaceBy = NULL) {
    templates <- gsub(pkgName, pattern, gsub("\\.rda$", ".Rd",
                                             list.files(file.path(pkgPath, pkgName, "data"))))

    for (i in templates) {
        if (is.null(replaceBy)) {
            rdName <- file.path(pkgPath, pkgName, "man",
                                gsub(pattern, pkgName, i))
        } else {
            rdName <- file.path(pkgPath, pkgName, "man",
                                gsub(pattern, replaceBy, i))
        }
        options(show.error.messages = FALSE)
        copied <- try(copySubstitute(file.path(.path.package("AnnBuilder"),
                                               "templates", i), rdName, repList, "#"))
        options(show.error.messages = TRUE)
        if (inherits(copied, "try-error")) {
            warning(paste("Can't copy", file.path(.path.package("AnnBuilder"),
                                                  "templates", i)))
        }
    }
}
