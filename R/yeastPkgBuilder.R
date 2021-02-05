# This function builds a data package for yeast geno.

yeastPkgBuilder <- function(pkgName, pkgPath, base="", version, author,
                            fromWeb=TRUE, lazyLoad=TRUE) {
    require("GO.db", quietly=TRUE) || stop("GO is needed to build data package")
    srcUrls = getSrcUrl("ALL", organism="Saccharomyces cerevisiae")
    oriUrls <- unlist(getOption("AnnBuilderPublicDataUrls"))
    
    srcObjs <- list()
    getUniCol4Yeast <- function(){
        if (base != "") {
            return(c("ORF", "CHR", "DESCRIPTION", "GENENAME"))
        } else {
            return(c("CHR", "DESCRIPTION", "GENENAME"))
        }
    }
    
    getMultCol4Yeast <- function() return(c("PMID", "ALIAS"))

    srcObjs[["YG"]] <- YG(srcUrl=srcUrls[["YG"]], fromWeb=fromWeb)
    if(!is.na(srcUrls[["KEGG"]])){
        srcObjs[["KEGG"]] <- KEGG(srcUrl = srcUrls[["KEGG"]],
                                  organism = "Saccharomyces cerevisiae",
                                  built = ifelse(fromWeb, getSrcBuilt("YG"),
                                  "Not available"), fromWeb = fromWeb,
                                  kegggenomeUrl = srcUrls[["KEGGGENOME"]])
    }
    if(!is.na(srcUrls[["GO"]])){
        srcObjs[["GO"]] <- GO(srcUrl = srcUrls[["GO"]],
                              built = ifelse(fromWeb, getSrcBuilt("GO"),
                              "Not Available"), fromWeb = fromWeb)
    }
    if(!is.na(srcUrls[["YEAST"]])){
        srcObjs[["YEAST"]] <- YEAST(srcUrl = srcUrls[["YEAST"]],
                              built = ifelse(fromWeb, getSrcBuilt("YEAST"),
                              "Not Available"), fromWeb = fromWeb)
    }
    ## to fix the "gopher5" error
    oriObjs <- srcObjs
    for(comp in names(srcObjs)){
        srcUrl(oriObjs[[comp]]) <- oriUrls[[toupper(comp)]]
    }

    makeSrcInfo()
    createEmptyDPkg(pkgName, pkgPath, force = TRUE)
    cat("package directory structure initialized\n")
    annotation <- yeastAnn(base)
    annotation <- as.matrix(annotation)
    colnames(annotation) <- toupper(colnames(annotation))
    cat("annotation data read complete\n")

    #goNCat <- mapGO2Category(goData)
    #annotation[, "GO"] <- as.vector(unlist(sapply(annotation[, "GO"],
    #                      nameGOByCat, goCat = goNCat), use.names = FALSE))
    # Write elements that has one to one mappings
    for(i in getUniCol4Yeast()){
        saveMat(annotation[, c("PROBE", i)], pkgName = pkgName,
                pkgPath = pkgPath, envName = i)
    }
    cat("1:1 PROBE mappings complete\n")
    # Write elements that has one to many mappings
    for(i in getMultCol4Yeast()){
        saveMat(annotation[, c("PROBE", i)], fun = splitEntry,
                pkgName = pkgName, pkgPath = pkgPath, envName = i)
    }
    cat("1:many PROBE mappings complete\n")

    
    ## YEASTCOMMON2SYSTEMATIC environment
    if (base == ""){
        revGN  <- annotation[, c("GENENAME","PROBE")]
        revAli <- annotation[, c("ALIAS","PROBE")]
        colnames(revGN) <- colnames(revAli) <- c("COMMON", "SYSTEMATIC")
        revAnno <- rbind(revGN, revAli)
        revAnno <- revAnno[!is.na(revAnno[, "COMMON"]), ]

        semiIndex <- grep(";", as.vector(revAnno[, "COMMON"]))
        listNonSys <- lapply(as.vector(revAnno[, "COMMON"]), function(x) unique(unlist(strsplit(x, ";"))))
        listNonSys <- lapply(listNonSys, function(x) ifelse(length(x)==0, return(NA), return(x)))

        eachRep <- sapply(listNonSys, length)
        revSYSTEMATIC <- rep(as.vector(revAnno[, "SYSTEMATIC"]), times=eachRep)
        revNonSys <- unlist(listNonSys)

        revAnno <- cbind(revNonSys, revSYSTEMATIC)
        colnames(revAnno) <- c("COMMON", "SYSTEMATIC")
        
        revAnnoSplit <- split(revAnno[, "SYSTEMATIC"], revAnno[, "COMMON"])
        revAnnoSplitCollapsed <- lapply(revAnnoSplit, function(x) paste(x[!is.na(x)], collapse=";"))
        revAnnoSplitCollapsed <- lapply(revAnnoSplitCollapsed, function(x) ifelse(x=="", return(NA), return(x)))
        revAnnoMat <- as.matrix(revAnnoSplitCollapsed)
        revAnnoMat <- cbind(rownames(revAnnoMat), revAnnoMat)
        colnames(revAnnoMat) <- c("COMMON", "SYSTEMATIC")

        saveMat(revAnnoMat, fun = splitEntry,
                pkgName=pkgName, pkgPath=pkgPath, envName="COMMON2SYSTEMATIC")
        cat("COMMON2SYSTEMATIC complete\n")
    }

    fun <- function(x) twoStepSplit(x, asNumeric = TRUE)
    saveMat(annotation[, c("PROBE", "CHRLOC")], fun = fun,
            pkgName = pkgName, pkgPath = pkgPath, envName = "CHRLOC")
    cat("CHRLOC complete\n")

    if(base==""){
        yeastDomain <- parseData(srcObjs[["YEAST"]])
        annotation <- merge(annotation, yeastDomain, by="PROBE", all.x=TRUE)
        annotation[annotation == ""] <- NA
        annotation <- as.matrix(annotation)

        ##fun <- function(x) twoStepSplit(x, asNumeric = FALSE)
        domainName <- setdiff(colnames(yeastDomain), "PROBE")
        for(domain in domainName){
            saveMat(annotation[, c("PROBE", domain)], fun = splitEntry,
                    pkgName = pkgName, pkgPath = pkgPath, envName = domain)
            cat(paste(domain, "complete\n"))
        }

        rejectORF(pkgPath=pkgPath)
        cat("rejectORF complete\n")
    }

    goCategories <- unlist(eapply(GOTERM, function(x) x@Ontology))
    goList <- getList4GO(goCategories,
                         sapply(annotation[, "GO"], twoStepSplit))
    names(goList) <- annotation[, "PROBE"]
    saveList(goList, pkgName, pkgPath, "GO")
    writeReverseMap(annotation[, c("PROBE", "PMID")], pkgName, pkgPath)

    if(!is.na(srcUrls[["KEGG"]])){
        if(base == ""){
            xkey <- "PROBE"
        }else{
            xkey <- "ORF"
        }
        # Adding KEGG pathway data to annoataion
        pathNEnzyme <- mapLL2ECNPName(srcObjs[["KEGG"]])
        orfNEC <- pathNEnzyme$llec
        orfNPath <- pathNEnzyme$llpathname
        colnames(orfNEC) <- c("ORF", "ENZYME")
        colnames(orfNPath) <- c("ORF", "PATH")
        annotation <- merge(annotation, orfNPath, by.x = xkey,
                            by.y = "ORF", all.x = TRUE)
        annotation <- merge(annotation, orfNEC, by.x = xkey,
                            by.y = "ORF", all.x = TRUE)
        annotation <- as.matrix(annotation)
        saveMat(annotation[, c("PROBE", "ENZYME")], fun = splitEntry,
            pkgName = pkgName, pkgPath = pkgPath, envName = "ENZYME")
        saveMat(annotation[, c("PROBE", "PATH")], fun = twoStepSplit,
            pkgName = pkgName, pkgPath = pkgPath, envName = "PATH")
        writeReverseMap(annotation[, c("PROBE", "PATH", "ENZYME")],
                        pkgName, pkgPath)
        cat("KEGG part complete\n")
    }
    go2Probe <- reverseMap4GO(annotation[, c("PROBE", "GO")], type = "ll2GO")
    colnames(go2Probe) <- c("GO", "GO2PROBE")
    go2ProbeEnv <- saveMat(go2Probe, fun=twoStepSplit,
                           pkgName=pkgName, pkgPath=pkgPath, envName="GO2PROBE")
    go2AllEnv <- createGo2AllProbesEnv(go2ProbeEnv)
    lockEnvironment(go2AllEnv, bindings=TRUE)
    envName <- paste(pkgName, "GO2ALLPROBES", sep="")
    assign(envName, go2AllEnv)
    fName <- file.path(pkgPath, pkgName, "data", paste(envName, ".rda", sep=""))
    save(list=envName, file=fName, compress=TRUE)
    if(base == ""){
        baseNew <- tempfile()
        write.table(getProbe2SGD(yGenoUrl=getSrcUrl("YG")), baseNew, sep = "\t", row.names = FALSE,
                    col.names = FALSE, quote = FALSE)
    }else{
        baseNew <- base
    }
    chrLengths <- getChrLengths("yeast")
    writeOrganism(pkgName, pkgPath, "yeast")
    writeChrLength(pkgName, pkgPath, chrLengths)
    ##repList <- getRepList("YEAST", srcObjs)
    repList <- getRepList("YEAST", oriObjs)
    repList[["PKGNAME"]] <- pkgName
    writeDocs(baseNew, pkgName, pkgPath, version, author,
              repList, "PKGNAME")

    if(base == ""){
        cat("Creating extra documents ...")
        for(doc in c("YEASTPFAM.Rd", "YEASTSMART.Rd", "YEASTINTERPRO.Rd")){
            copySubstitute(file.path(.path.package("AnnBuilder"),
                                     "templates", doc),
                           file.path(pkgPath, pkgName, "man", doc),
                           list(YEASTDOMAIN="Saccharomyces Genome Database",
                                DATE=date(),
                                URL=paste(srcUrl(srcObjs[["YEAST"]]), baseFile(srcObjs[["YEAST"]]), sep="")),
                           "#")
        }
        cat(" Done\n")
    }
    
    if (lazyLoad)
      makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
}
