athPkgBuilder <- function(
                          baseName=NULL,
                          pkgName,
                          pkgPath,
                          fileExt = list(
			                base = "Microarrays/Affymetrix/affy_ATH1_array_elements-2006-07-14.txt",
                            estAssign = "Genes/est_mapping/est.Assignment.Locus",
                            seqGenes = "Genes/TAIR_sequenced_genes",
                            go = "Ontologies/Gene_Ontology/ATH_GO_GOSLIM.20050827.txt",
                            aliases = "Genes/gene_aliases.20041105",
                            aracyc = "Pathways/aracyc_dump_20050412",
                            kegg = "/ath/ath_gene_map.tab",
                            pmid = "User_Requests/LocusPublished.08012006.txt"),
                          ncols = list(
                            base = 9,
                            estAssign = 7,
                            seqGenes = 4,
                            go = 12,
                            aliases = 4,
                            aracyc = 4,
                            kegg = 2,
                            pmid = 4),
                          cols2Keep = list(
                            base = c(1, 5),
                            estAssign = c(3, 6, 7),
                            seqGenes = c(1, 3, 4),
                            go = c(1, 5, 9),
                            aliases = c(1, 2),
                            aracyc = c(1, 3, 4),
                            kegg = c(1,2),
                            pmid = c(1, 4)),
                          colNames = list(
                            base = c("PROBE", "ACCNUM"),
                            estAssign = c("CHRLOC", "ORI",  "ACCNUM"),
                            seqGenes = c("ACCNUM", "CHR", "GENENAME"),
                            go = c("ACCNUM", "GO", "EVID"),
                            aliases = c("ACCNUM", "SYMBOL"),
                            aracyc = c("ARACYC", "ENZYME", "ACCNUM"),
                            kegg = c("ACCNUM", "PATH"),
                            pmid = c("ACCNUM", "PMID")),
                          indexby = "PROBE",
                          version,
                          author,
                          lazyLoad=TRUE) {
    require("GO.db", quietly=TRUE) || stop("GO is needed to build data package")
    baseUrl <- getSrcUrl("AT")
    #oriUrls <- unlist(getOption("AnnBuilderPublicDataUrls"))
    
    #srcObjs <- getSrcObjs(sapply("GO", getSrcUrl,
    #                      organism = "Arabidopsis thaliana"),
    #                      baseName = "", "", baseMapType = "gb")
    #oriObjs <- srcObjs
    #for(comp in names(srcObjs)){
    #    srcUrl(oriObjs[[comp]]) <- oriUrls[[toupper(comp)]]
    #}
    if (is.null(baseName)) {
        base <- readAthData(baseUrl, fileExt[["base"]], cols2Keep[["base"]],
                            colNames[["base"]], ncols[["base"]])
        base[,"ACCNUM"] <- sub("^no_match$", "NA", base[,"ACCNUM"])
        base <- base[-1,]
	if (indexby=="PROBE") {
	    baseName <- tempfile()
  	    write.table(base, baseName, sep = "\t", col.names = FALSE, row.names = FALSE, quote=F)
	}
    } else {        
        base <- read.table(baseName, sep = "\t", header = FALSE, as.is = TRUE, comment.char="", stringsAsFactors=FALSE)
        base[,2] <- toupper(base[,2])
        base <- mergeDupMatByFirstCol(base)
        colnames(base) <- c("PROBE", "ACCNUM")
    }
    multiIndex <- grep(";", base[,"ACCNUM"])
    # tempBase contains all the AGI locus that we need to annotate
    if (indexby=="PROBE") {
        tempBase <- base[-multiIndex, ]
    } else if(indexby=="ACCNUM") {
        tempBase <- unique(unlist(strsplit(base[,"ACCNUM"], ";")))
        tempBase <- tempBase[tempBase!="NA"]
        tempBase <- as.data.frame(tempBase)
        colnames(tempBase) <- "ACCNUM" 
    } else {
        stop("Invalid 'indexby', must be either 'PROBE' or 'ACCNUM'.")
    }
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
    seqGenes <- readAthData(baseUrl, fileExt[["seqGenes"]],
                          cols2Keep[["seqGenes"]], colNames[["seqGenes"]],
                         ncols[["seqGenes"]])
    chrom <- seqGenes[2:nrow(seqGenes), 1:2] #acc and chrom
    chrom <- getOneMap(chrom, "ACCNUM")
    athData <- merge(athData, chrom, by = "ACCNUM", all.x = TRUE)
    # Get descriptions for genes
    genename <- seqGenes[2:nrow(seqGenes), c(1,3)] #acc and description
    colnames(genename) <- c("ACCNUM", "GENENAME")
    genename <-  mergeDupMatByFirstCol(genename)
    athData <- merge(athData, genename, by="ACCNUM", all.x = TRUE)
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
    #Get AraCyc data
    aracyc <- readAthData(baseUrl, fileExt[["aracyc"]],
                          cols2Keep[["aracyc"]], colNames[["aracyc"]],
                            ncols[["aracyc"]])
    aracyc[, "ACCNUM"] <- toupper(aracyc[, "ACCNUM"])
    subEnzyme = matrix(lapply(aracyc[,2], function(x) strsplit(x,"\\|\\|")[[1]][1]), ncol=1)
    aracyc = data.frame(as.character(aracyc[,1]), as.character(subEnzyme), as.character(aracyc[,3]))
    colnames(aracyc) <- c("ARACYC", "ENZYME", "ACCNUM")

    aracyc <- mergeDupMatByFirstCol(aracyc[, c("ACCNUM",
              colnames(aracyc)[colnames(aracyc) != "ACCNUM"])])
    athData <- merge(athData, aracyc, by = "ACCNUM", all.x = TRUE)
    
    #Get KEGG data
    kegg <- readAthData(getSrcUrl("KEGG"), fileExt[["kegg"]],
                          cols2Keep[["kegg"]], colNames[["kegg"]],
                            ncols[["kegg"]])
    kegg[, "ACCNUM"] <- toupper(kegg[, "ACCNUM"])
    #kegg <- cbind(kegg[,"ACCNUM"], strsplit(kegg[, "PATH"], " "))
    kegg[,"PATH"] <- gsub(" ", ";", kegg[,"PATH"])
    colnames(kegg) <- c("ACCNUM", "PATH")
    athData <- merge(athData, kegg, by = "ACCNUM", all.x = TRUE)
    # Get PMID data
    pmid <- readAthData(baseUrl, fileExt[["pmid"]], cols2Keep[["pmid"]],
                        colNames[["pmid"]], ncols[["pmid"]])
    naPM <- which(pmid[,"PMID"]=="")
    if(length(naPM)>0)
	pmid<- pmid[-naPM,]
    pmid <- mergeDupMatByFirstCol(pmid[, c("ACCNUM", "PMID")])
    athData <- merge(athData, pmid, by = "ACCNUM", all.x = TRUE)
    
    # deal with the case when one probeset hits multiple AGI locus
    if (indexby == "PROBE") {
        athData <- cbind(athData, NA)
        colnames(athData)[12] <- "MULTIHIT"
        base <- base[multiIndex,]
        base <- cbind("multiple", base[,"PROBE"], "multiple", "multiple", "multiple",
            "multiple", "multiple", "multiple", "multiple","multiple", "multiple",
            base[,"ACCNUM"])
        colnames(base) <- colnames(athData)
        # assume athData has 12 columns: ACCNUM, PROBE, CHRLOC, CHR, 
        # GENENAME, GO, SYMBOL, ARACYC, ENZYME, PATH, PMID, MULTIHIT
        athData <- rbind(athData, base)
    } else { # indexby="ACCNUM"
        colnames(athData)[1] <- "PROBE"
    } 
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
    if(indexby=="ACCNUM") {
        accnumList <- strsplit(base[,"ACCNUM"], ";")
	base <- cbind(base[,"PROBE"], accnumList)
        saveMat(base, pkgName=pkgName, pkgPath=pkgPath, envName="ACCNUM")
    }
    ## to fix the "gopher5" problem due to using local missor
    oriUrls <- unlist(getOption("AnnBuilderPublicDataUrls"))
    baseUrl <- oriUrls[["AT"]]
    repList <- list(LLSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("base", "seqGenes", "pmid")], "}",
                        sep = "", collapse = "\n"), "\nBuilt: N/A",
                        sep = ""),
                    GOSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("go")], "}", sep = "",
                        collapse = "\n"), "\nBuilt: N/A", sep = ""),
                    GPSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("estAssign")], "}", sep = "",
                        collapse = "\n"), "\nBuilt: N/A", sep = ""),
                    KEGGSOURCE = paste(paste("Tair: \\\\url{", baseUrl,
                        fileExt[c("kegg")], "}", sep = "",
                        collapse = "\n"), "\nBuilt: N/A", sep = ""),
                    DATE = date())
     repList[["PKGNAME"]] <- pkgName
     options(show.error.messages = FALSE)
     organism <- "Arabidopsis thaliana"
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
        temp <- temp[grep("^AT", temp)]
        if(length(temp) == 0){
            return(NA)
        }else{
            return(temp)
        }
    }
    for(i in 1:nrow(pmid)){
        ##print(pmid[i,])
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
                          sep = "\t", as.is = TRUE, fill=TRUE, stringsAsFactors=FALSE,
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

getFileExt <- function(chipName="ATH1", verbose=FALSE) {
    if(!(chipName %in% c("ATH1", "AG")))
        stop("chipName should be either 'ATH1' or 'AG'.")
    require("RCurl", quietly=TRUE)||stop("RCurl is needed to create FileExt.")
    baseUrl <- getSrcUrl("AT")
    curl <- getCurlHandle()
    # get url of gene-to-locus mapping
    locusUrl <- "Microarrays/Affymetrix/"
    locusFile <- paste("^.*(affy_", chipName, "_array_elements-[^.]*\\.txt).*$", sep="")
    locusFile <- sub(locusFile, "\\1", getURL(paste(baseUrl, locusUrl, sep=""), ftp.use.epsv=0, verbose=verbose))
    locusUrl <- paste(locusUrl, locusFile, sep="")
    # get url of go
    goUrl <- "Ontologies/Gene_Ontology/"
    goFile <- sub("^.*(ATH_GO_GOSLIM\\.[0-9]*\\.txt).*$", "\\1",
                getURL(paste(baseUrl, goUrl, sep=""), ftp.use.epsv=0, verbose=verbose))
    goUrl <- paste(goUrl, goFile, sep="")
    # get url of aliases
    aliasesUrl <- "Genes/"
    aliasesFile <- sub("^.*(gene_aliases\\.[0-9]*).*$", "\\1",
                    getURL(paste(baseUrl, aliasesUrl, sep=""), ftp.use.epsv=0, verbose=verbose))
    aliasesUrl <- paste(aliasesUrl, aliasesFile, sep="")
    # get url of aracyc
    aracycUrl <- "Pathways/"
    aracycFile <- sub("^.*(aracyc_dump_[0-9]*).*$", "\\1",
                    getURL(paste(baseUrl, aracycUrl, sep=""), ftp.use.epsv=0, verbose=verbose))
    aracycUrl <- paste(aracycUrl, aracycFile, sep="")
    # get url of pmid
    pmidUrl <- "User_Requests/"
    pmidFile <- sub("^.*(LocusPublished\\.[0-9]*(\\.txt)?).*$", "\\1",
                    getURL(paste(baseUrl, pmidUrl, sep=""), ftp.use.epsv=0, verbose=verbose))
    pmidUrl <- paste(pmidUrl, pmidFile, sep="")
    # assemble fileExt
    fileExt <- list(
        base = locusUrl,
        estAssign = "Genes/est_mapping/est.Assignment.Locus",
        seqGenes = "Genes/TAIR_sequenced_genes",
        go = goUrl,
        aliases = aliasesUrl,
        aracyc = aracycUrl,
        kegg = "/ath/ath_gene_map.tab",
        pmid = pmidUrl
    )
    fileExt
}
        
    
     
    

