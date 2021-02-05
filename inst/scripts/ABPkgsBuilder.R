ABPkgsBuilder <- function(basePath ="/scratch/homes/jzhang/affymappings",
                          baseMapType = "gbNRef",
                          pkgPath = "/scratch/homes/jzhang/1.6.2",
                          version = "1.6.2",
                          author = list(author = "Jianhua Zhang",
                          maintainer = "jzhang@jimmy.harvard.edu"),
                          pattern = "_.*\\.txt$", makeXML = TRUE){


    makeSrcInfo()
    srcs <- getPkgFileNames(basePath, pattern = pattern)
    orgs <- unique(unlist(sapply(names(srcs), getOrganism),
                          use.names = FALSE))
    srcObjs <- getSrcObjs(sapply(c("LL", "GO"), getSrcUrl,
                          organism = ""), baseName = "", "", baseMapType)
    goData <- readData(srcObjs[["go"]], xml = TRUE)
    for(j in orgs){
        organism <- j
        srcObjs[c("ug", "kegg", "gp")] <-
            getSrcObjs(sapply(c("UG", "KEGG", "GP"), getSrcUrl,
                       organism = organism), baseName = "",
                       organism, baseMapType)

        # Get strand data by processing data from GoldenPath
        options(show.error.messages = FALSE)
        strand <- try(getStrand(srcObjs[["gp"]]))
        options(show.error.messages = TRUE)
        if(inherits(strand, "try-error")){
            warning(paste("Failed to parse Golden Path data because of:\n\n",
                       strand))
            # If failed, create an empty matrix for strand so that the
            # subsequent steps can go through
            strand <- matrix(nrow = 0, ncol = 2)
            colnames(strand) <- c("LOCUSID", "CHRLOC")
        }
        # Get KEGG pathway data
        options(show.error.messages = FALSE)
        pathNEnzyme <- try(mapLL2ECNPName(srcObjs[["kegg"]]))
        options(show.error.messages = TRUE)
        if(inherits(pathNEnzyme, "try-error")){
            warning(paste("Failed to parse KEGG data because of:\n\n",
                       pathNEnzyme))
            path <-  matrix(NA, ncol = 2)
            colnames(path) <- c("LOCUSID", "PATH")
            ec <-  matrix(NA, ncol = 2)
            colnames(ec) <- c("LOCUSID", "ENZYME")
            pathNEnzyme <- list(llpathname = path, llec = ec)
        }

        for(i in names(srcs)){
            print(paste("Processing package", i))
            baseFile(srcObjs[["ll"]]) <- srcs[[i]][["ID"]]
            baseFile(srcObjs[["ug"]]) <- srcs[[i]][["ID"]]
            parser(srcObjs[["ll"]]) <- getBaseParsers("gb")["LL"]
            otherSrc <- srcs[[i]][names(srcs[[i]]) != "ID"]
            pkgName <- i
            baseName <-  srcs[[i]][["ID"]]
            unified <- getUniMappings(srcs[[i]][["ID"]], srcObjs[["ll"]],
                                      srcObjs[["ug"]], otherSrc, baseMapType)
            if(baseMapType == "ll"){
                parser(srcObjs[["ll"]]) <- file.path(.path.package(
                                     "AnnBuilder"), "scripts", "llParser4LL")
            }else{
                parser(srcObjs[["ll"]]) <- file.path(.path.package(
                                     "AnnBuilder"), "scripts", "llParser")
            }
            baseFile(srcObjs[["ll"]]) <- unified
            annotation <- getAnnData(srcObjs["ll"])
            annotation <- merge(annotation, strand, by = "LOCUSID",
                                all.x = TRUE)
            annotation <- merge(annotation, pathNEnzyme$llpathname,
                                by = "LOCUSID", all.x = TRUE)
            annotation <- merge(annotation, pathNEnzyme$llec,
                                by = "LOCUSID", all.x = TRUE)
            annotation <- as.matrix(annotation)
            annotation[annotation == "" | annotation == "NA"] <- NA

            if (makeXML) {
                multC <- colnames(annotation)[is.element(colnames(annotation),
                                                        getMultiColNames())]
                typeC <- colnames(annotation)[is.element(colnames(annotation),
                                                        getTypeColNames())]
                XMLOut <- file.path(pkgPath,
                                    paste(pkgName, ".xml", sep = ""))
                fileToXML(targetName = pkgName, outName = XMLOut,
                          inName = annotation, colNames = "",
                          idColName = "PROBE",
                          multColNames = multC, typeColNames = typeC,
                          isFile = FALSE, version = version)
            }

            createEmptyDPkg(pkgName = pkgName, pkgPath, force = TRUE)
            writeAnnData2Pkg(annotation, pkgName, pkgPath)
            revNames <- intersect(colnames(annotation),
                                  c("PMID", "PATH", "ENZYME"))
            if(length(revNames) != 0){
                writeReverseMap(annotation[, c("PROBE", revNames)],
                                pkgName, pkgPath)
            }
            if(!all(is.na(annotation[, "GO"]))){
                goList <- getList4GO(getNamedCat(goData),
                                 sapply(annotation[, "GO"], twoStepSplit))
                names(goList) <- annotation[, "PROBE"]
                saveList(goList, pkgName, pkgPath, "GO")
                go2Probe <- mapGO2Probe(srcObjs[["ll"]], baseMapType)
                saveMat(go2Probe, fun = twoStepSplit, pkgName = pkgName,
                    pkgPath = pkgPath, envName = "GO2PROBE")
                go2ALL <- mapGO2AllProbe(go2Probe, goData, "",
                                     sep = ";", all = TRUE)
                saveMat(cbind(names(go2ALL), go2ALL), fun = twoStepSplit,
                    pkgName = pkgName, pkgPath = pkgPath,
                    envName = "GO2ALLPROBES")

                if (makeXML) {
                    go2ALL <- cbind(names(go2ALL), go2ALL)
                    colnames(go2ALL) <- c("GO", "GO2ALLPROBES")
                    if (!is.null(go2ALL) && !is.null(go2Probe)) {
                        mergedGO <- as.matrix(merge(as.matrix(go2ALL),
                              as.matrix(go2Probe), by = "GO", all.x = TRUE))
                        XMLByNum <- file.path(pkgPath, paste(pkgName,
                                                     "ByNum.xml", sep = ""))
                        fileToXML(targetName = "GOByNum", outName = XMLByNum,
                                  inName = mergedGO, colNames = "",
                                  idColName = "GO",
                               multColNames = c("GO2ALLPROBES", "GO2PROBE"),
                                  typeColNames = "",
                                  isFile = FALSE, version = version)
                    }
                }
            }

            repList <- getRepList("all", resumeSrcUrl(srcObjs))
            repList[["PKGNAME"]] <- pkgName
            options(show.error.messages = FALSE)
            chrLengths <- try(getChrLengths(organism))
            options(show.error.messages = TRUE)
            if(!inherits(chrLengths, "try-error")){
                writeChrLength(pkgName, pkgPath, chrLengths)
            }
            writeOrganism(pkgName, pkgPath, organism)
            writeDocs(baseName, pkgName, pkgPath, version, author,
                      repList, "PKGNAME")
            makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
        }
        cleanSrcObjs(srcObjs["ug"])
    }
    cleanSrcObjs(srcObjs["ll"])
}

getPkgFileNames <- function(basePath, pattern = "_.*\\.txt$"){
    # Create a list for holding the data
    idNSrcs <- list()

    pkgNames <- list.files(basePath, pattern = pattern, full.names = TRUE)

    for(i in unique(gsub("_.*\\.txt", "", pkgNames))){
        files4OnePkg <- pkgNames[grep(paste("^", i, "_.*", sep = ""),
                                     pkgNames)]
        if(length(grep("_ID", files4OnePkg)) > 0){
            files4OnePkg <- c(files4OnePkg[grep("_ID", files4OnePkg)],
                              files4OnePkg[-grep("_ID", files4OnePkg)])
            names(files4OnePkg) <- c("ID",
                                  LETTERS[1:length(files4OnePkg) - 1])
            idNSrcs[[basename(i)]] <- files4OnePkg
        }
    }
    return(idNSrcs)
}

getOrganism <- function(orgChar){

    orgChar <- substr(orgChar, 1, 2)

    switch(toupper(orgChar),
           HG = ,
           HU = return("Homo sapiens"),
           MO = ,
           MG = return("Mus musculus"),
           RA = ,
           RG = return("Rattus norvegicus"),
           AT = return("Arabidopsis thaliana"),
           DR = return("Drosophila melanogaster"),
           ZE = return("Danio rerio"),
           XE = return("Xenopus laevis"),
           CE = return("Caenorhabditis elegans"),
           stop(paste("Unknown entry for orgChar:", orgChar)))
}

getPkgNameByOrg <- function(pkgPath, organism){
    switch(toupper(organism),
           HUMAN = init <- "h",
           MOUSE = init <- "mg\|mo",
           RAT = init <- "rg\|ra",
           ARABIDOPSIS = init <- "at",
           FLY = init <- "dros",
           stop(paste("Ubknown organism", organism)))

    pkgNames <- list.files(pkgPath, pattern = init)

    return(setdiff(pkgNames, pkgNames[grep("\\.xml$", pkgNames)]))
}

putDPkgsToWeb <- function(pkgPath = "/scratch/homes/jzhang/datapkgs",
                    webPath = "/scratch/homes/jzhang/bioC/webstuff/data",
                    organism = "human", version = "1.2.1", what = "gz"){

    if(organism == "GO" || organism == "KEGG"){
        destDir <- file.path(webPath, "misc", "files", version)
    }else{
        destDir <- file.path(webPath, organism, "files", version)
    }

    if(!file.exists(destDir)){
        dir.create(destDir)
    }

    system(paste("cp", file.path(pkgPath, paste(substr(organism, 1, 1),
                           paste("*.", what, sep = ""), sep = "")), destDir))
}

writeDataHTML <- function(cin =
    "/scratch/homes/jzhang/madman/Rpacks/AnnBuilder/inst/scripts/data.html",
    cout = "/scratch/homes/jzhang/bioC/webstuff/data/data.html",
                          symbolValues = list(VERSION = "1.2.1")){
    copySubstitute(cin, cout, symbolValues)
}

checkNBuildPkgs <- function(pkgPath, check = TRUE, exclude = NULL){
    curDir <- getwd()
    setwd(pkgPath)
    pkgs <- list.files(pkgPath)
    if(!is.null(exclude)){
        pkgs <- setdiff(pkgs, exclude)
    }
    # Get rid of any file with an extension
    pkgs <- pkgs[-grep("\\.", pkgs)]
    if(check){
        for(i in pkgs){
            system(paste("R CMD INSTALL", i))
            library(i, character.only = TRUE)
            # Check for probe missmatch for probe based envs
            pBased <- getPBasedEnv()
            for(k in pBased[-1]){
                if(any(!is.element(ls(get(paste(i, pBased[1], sep = ""))),
                                   ls(get(paste(i, k, sep = "")))))){
                    stop(paste("Checking for probe missmatch failed for",
                               pBased[1], "and", k))
                }
            }
            # checking for probe missmatch for non-probe based envs
            for(k in getNPEnvs()){
                probes <- ls(get(paste(i, "ACCNUM", sep = "")))
                toMatch <- unique(unlist(multiget(ls(get(paste(i, k,
                                sep = ""))), get(paste(i, k, sep = ""))),
                                         use.names = FALSE))
                toMatch <- toMatch[!is.na(toMatch)]
                if(any(!is.element(toMatch, probes))){
                    stop(paste("Checking for", paste(i, k, sep = ""),
                               "failed"))
                }
            }
        }
    }
    # If checking is successful, build the packages
    for(i in pkgs){
        system(paste("R CMD build -f", i))
    }
    setwd(curDir)
}

getNPEnvs <- function(){
    return(c("ENZYME2PROBE", "GO2ALLPROBES", "GO2PROBE", "PATH2PROBE",
           "PMID2PROBE"))
}

getPBasedEnv <- function(){
    return(c("ACCNUM", "CHR", "CHRLOC" , "GENENAME", "GO", "HGID",
             "LOCUSID", "MAP", "PATH", "PMID", "SUMFUNC",
             "SYMBOL", "UNIGENE", "OMIM", "NP", "NM", "REFSEQ"))
}

getAllSrcData <- function(srcUrls, organisms){
    llSrc <- loadFromUrl(getSrcUrl("LL"))
    ugSrc <- list()
    gpSrc <- list()
    keggSrc <- list()
    hgid <- list()
    for(i in organisms){
        ugSrc[[i]] <- loadFromUrl(getSrcUrl("UG", i))
        gp <- GP(srcUrl = getSrcUrl("GP", organism = i), organism = i)
        options(show.error.messages = FALSE)
        strand <- try(getStrand(gp))
        options(show.error.messages = TRUE)
        if(inherits(strand, "try-error")){
            stop(paste("Failed to parse Golden Path data because of:\n\n",
                       strand))
        }else{
            gpSrc[[i]] <- strand
        }
        kegg <- KEGG(srcUrl = getSrcUrl("KEGG", organism = i), organism = i)
        keggSrc[[i]] <- mapLL2ECNPName(kegg)
        hgid[[i]] <- getLL2IntID(procHomoData(getSrcUrl("HG")), i)
    }
}
