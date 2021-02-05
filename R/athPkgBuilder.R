athPkgBuilder <- function(baseName, pkgName, pkgPath,
                   fileExt = list(estAssign =
                   "Genes/est_mapping/est.Assignment.Locus",
                   seqGenes = "Genes/TAIR_sequenced_genes",
                   go = "Ontologies/Gene_Ontology/ATH_GO_GOSLIM.20050507.txt",
                   aliases = "Genes/gene_aliases.20041105",
                   pathway = "Pathways/aracyc_dump_20050412",
                   pmid = "Genes/Gene_Anatomy/ATH_Anatomy.20040209.txt"),
                   ncols = list(estAssign = 7, seqGenes = 4, go = 12,
                   aliases = 4, pathway = 4, pmid = 15),
    cols2Keep = list(estAssign = c(3, 6, 7), seqGenes = c(1, 3),
                   go = c(1, 5, 9), aliases = c(1, 2),
                   pathway = c(1, 3, 4), pmid = c(6, 11)),
    colNames = list(estAssign = c("CHRLOC", "ORI",  "ACCNUM"),
                   seqGenes = c("ACCNUM", "CHR"),
                   go = c("ACCNUM", "GO", "EVID"),
                   aliases = c("ACCNUM", "SYMBOL"),
                   pathway = c("PATH", "ENZYME", "ACCNUM"),
                   pmid = c("PMID", "ACCNUM")),
    version,
    author,
    baseUrl="ftp://ftp.arabidopsis.org/home/tair/", lazyLoad=TRUE) {
    require("GO", quietly=TRUE) || stop("GO is needed to build data package")

    srcObjs <- getSrcObjs(sapply("GO", getSrcUrl,
                          organism = "Arabidopsis thaliana"),
                          baseName = "", "", baseMapType = "gb")
    base <- read.table(baseName, sep = "\t", header = FALSE, as.is = TRUE)
    tempBase <- cbind(base[,1], toupper(base[,2]))
    colnames(tempBase) <- c("PROBE", "ACCNUM")
    # Get the location and orientation data for genes
    athData <- readAthData(baseUrl, fileExt[["estAssign"]],
                           cols2Keep[["estAssign"]], colNames[["estAssign"]],
                           ncols[["estAssign"]])
    athData <- athData[2:nrow(athData),]
    athData <- athData[athData[,3] != "",]
    athData[, 2] <- ifelse(athData[,2] == "forward", "+", "-")
    athData <- rbind(getOneMap(athData[athData[, "ORI"] == "+", ], "ACCNUM"),
                     getOneMap(athData[athData[, "ORI"] == "-", ], "ACCNUM"))
    athData <- cbind(athData[,3],
               paste(athData[,2], athData[, 1], sep = ""))
    colnames(athData) <- c("ACCNUM", "CHRLOC")
    athData <-  mergeDupMatByFirstCol(athData)
    athData <- merge(tempBase, athData, by = "ACCNUM", all.x = TRUE)
    # Get chromosome number for genes
    chrom <- readAthData(baseUrl, fileExt[["seqGenes"]],
                          cols2Keep[["seqGenes"]], colNames[["seqGenes"]],
                         ncols[["seqGenes"]])
    chrom <- chrom[2:nrow(chrom),]
    chrom <- getOneMap(chrom, "ACCNUM")
    athData <- merge(athData, chrom, by = "ACCNUM", all.x = TRUE)
    # Get GO data
    go <- readAthData(baseUrl, fileExt[["go"]], cols2Keep[["go"]],
                      colNames[["go"]], ncols[["go"]])
    go <- cbind(go[,1], paste(go[, 2], go[, 3], sep = "@"))
    colnames(go) <- c("ACCNUM", "GO")
    go <- mergeDupMatByFirstCol(go)
    athData <- merge(athData, go, by = "ACCNUM", all.x = TRUE)
    # Get gene symbol
    symbol <- readAthData(baseUrl, fileExt[["aliases"]],
                          cols2Keep[["aliases"]], colNames[["aliases"]],
                          ncols[["aliases"]])
    
    colnames(symbol) <- c("ACCNUM", "SYMBOL")
    symbol <- mergeDupMatByFirstCol(symbol)
    athData <- merge(athData, symbol, by = "ACCNUM", all.x = TRUE)
    #Get pathway data
    pathways <- readAthData(baseUrl, fileExt[["pathway"]],
                          cols2Keep[["pathway"]], colNames[["pathway"]],
                            ncols[["pathway"]])
    pathways[, "ACCNUM"] <- toupper(pathways[, "ACCNUM"])
    subEnzyme = matrix(lapply(pathways[,2], function(x) strsplit(x,"\\|\\|")[[1]][1]), ncol=1)
    pathways = data.frame(as.character(pathways[,1]), as.character(subEnzyme), as.character(pathways[,3]))
    colnames(pathways) <- c("PATH", "ENZYME", "ACCNUM")

    pathways <- mergeDupMatByFirstCol(pathways[, c("ACCNUM",
              colnames(pathways)[colnames(pathways) != "ACCNUM"])])
    athData <- merge(athData, pathways, by = "ACCNUM", all.x = TRUE)

    # Get PMID data
    pmid <- readAthData(baseUrl, fileExt[["pmid"]], cols2Keep[["pmid"]],
                        colNames[["pmid"]], ncols[["pmid"]])
    pmid <- procPMIDData(pmid)
    pmid <- mergeDupMatByFirstCol(pmid[, c("ACCNUM", "PMID")])
    athData <- merge(athData, pmid, by = "ACCNUM", all.x = TRUE)
    athData <- as.matrix(athData)
    # Write data
    makeSrcInfo()
    createEmptyDPkg(pkgName = pkgName, pkgPath, force = TRUE)
    writeAnnData2Pkg(athData, pkgName, pkgPath)
    revNames <- intersect(colnames(athData), c("PMID", "PATH", "ENZYME"))
    if(length(revNames) != 0){
        writeReverseMap(athData[, c("PROBE", revNames)], pkgName, pkgPath)
    }
    if(!all(is.na(athData[, "GO"]))){
        goCategories <- unlist(eapply(GOTERM, function(x) x@Ontology))
        goList <- getList4GO(goCategories,
                             sapply(athData[, "GO"], twoStepSplit))
        names(goList) <- athData[, "PROBE"]
        saveList(goList, pkgName, pkgPath, "GO")
        go2Probe <- reverseMap4GO(athData[, c("PROBE", "GO")], type = "ll2GO")
        colnames(go2Probe) <- c("GO","PROBE")
        go2ProveEnv <- saveMat(go2Probe, fun = twoStepSplit, pkgName = pkgName,
                               pkgPath = pkgPath, envName = "GO2PROBE")
        go2AllEnv <- createGo2AllProbesEnv(go2ProveEnv)
        lockEnvironment(go2AllEnv, bindings=TRUE)
        envName <- paste(pkgName, "GO2ALLPROBES", sep="")
        assign(envName, go2AllEnv)
        fName <- file.path(pkgPath, pkgName, "data", paste(envName, ".rda", sep=""))
        save(list=envName, file=fName, compress=TRUE)
    }
    repList <- list(LLSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("aliases", "seqGenes", "pmid")], "}",
                        sep = "", collapse = "\n"), "\nBuilt: N/A",
                        sep = ""),
                    GOSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("go")], "}", sep = "",
                        collapse = "\n"), "\nBuilt: N/A", sep = ""),
                    GPSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("estAssign")], "}", sep = "",
                        collapse = "\n"), "\nBuilt: N/A", sep = ""),
                    KEGGSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("pathway")], "}", sep = "",
                        collapse = "\n"), "\nBuilt: N/A", sep = ""),,
                    DATE = date())
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
    if (lazyLoad)
      makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)

     return(invisible())
}

getOneMap <- function(map, keyCol){
    return(map[!duplicated(map[, keyCol]),])
}

procPMIDData <- function(pmid){
    processed <- matrix(ncol = 2, nrow = 0)
    procACC <- function(acc){
        temp <- unlist(strsplit(acc, "\\|"))
        temp <- temp[grep("AT", temp)]
        if(length(temp) == 0){
            return(NA)
        }else{
            return(temp)
        }
    }
    for(i in 1:nrow(pmid)){
        print(pmid[i,])
        processed <- rbind(processed, cbind(gsub(".*PMID:(.*$)","\\1",
                                           pmid[i,1]), procACC(pmid[i,2])))
    }
    colnames(processed) <- colnames(pmid)
    return(processed)
}

getSrcObjs4Ath <- function(){
    srcObjs <- list()
    srcObjs <- getSrcObjs(sapply("GO", getSrcUrl,
                          organism = "Arabidopsis thaliana"),
                          baseName = "", "", baseMapType = "gb")
    return(srcObjs)
}

readAthData <- function(baseUrl, ext, col2Keep, colNames, ncols){
    options(show.error.messages = FALSE)
    targetURL <- loadFromUrl(paste(baseUrl, ext, sep = ""))
    #athData <- try(matrix(scan(targetURL, what = "", sep = "\t", quote = "",
    #                    strip.white = TRUE, quiet = TRUE, fill = TRUE),
    #                      ncol = ncols, byrow = TRUE))
    athData <- try(read.table(paste(baseUrl, ext, sep = ""),
                          sep = "\t", as.is = TRUE, fill=TRUE,
                        quote = "", comment.char = ""))
    options(show.error.messages = TRUE)
    if(inherits(athData, "try-error")){
        stop(paste("Failed to obtain data from",
                   paste(baseUrl, ext, sep = "")))
    }else{
        athData <- athData[, col2Keep]
        if(!missing(colNames)){
            colnames(athData) <- colNames
        }
        return(athData)
    }
}

mergeDupMatByFirstCol <- function(dupMat, sep = ";"){
    # Retain column names
    colNames <- colnames(dupMat)

    merged <- t(sapply(split.data.frame(dupMat, factor(dupMat[,1])),
                   function(x){
                       apply(x, 2, FUN = function(y){
                           return(paste(unique(y), sep = "", collapse = sep))
                       })
                   }))
   if(!is.null(colNames)){
       colnames(merged) <- colNames
   }
   return(merged)
}

