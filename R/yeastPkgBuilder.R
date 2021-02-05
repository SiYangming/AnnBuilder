# This function builds a data package for yeast geno.

yeastPkgBuilder <- function(pkgName, pkgPath, base="", version, author,
                            fromWeb=TRUE, lazyLoad=TRUE) {
    require("GO", quietly=TRUE) || stop("GO is needed to build data package")
    srcUrls = getSrcUrl("ALL", organism="Saccharomyces cerevisiae")
    
    srcObjs <- list()
    
    getUniCol4Yeast <- function(){
        if (base != "") {
            return(c("PMID", "ORF", "CHR", "DESCRIPTION", "GENENAME"))
        } else {
            return(c("PMID", "CHR", "DESCRIPTION", "GENENAME"))
        }
    }
    
    getMultCol4Yeast <- function() return("ALIAS")

    srcObjs[["YG"]] <- YG(srcUrl=srcUrls[["YG"]], fromWeb=fromWeb)
    if(!is.na(srcUrls[["KEGG"]])){
        srcObjs[["KEGG"]] <- KEGG(srcUrl = srcUrls[["KEGG"]],
                                  organism = "Saccharomyces cerevisiae",
                                  built = ifelse(fromWeb, getSrcBuilt("YG"),
                                  "Not available"), fromWeb = fromWeb)
    }
    if(!is.na(srcUrls[["GO"]])){
        srcObjs[["GO"]] <- GO(srcUrl = srcUrls[["GO"]],
                              built = ifelse(fromWeb, getSrcBuilt("GO"),
                              "Not Available"), fromWeb = fromWeb)
    }

    makeSrcInfo()
    createEmptyDPkg(pkgName, pkgPath, force = TRUE)
    cat("package directory structured initialized\n")
    annotation <- yeastAnn(base)
    annotation <- as.matrix(annotation)
    annotation[annotation == ""] <- NA
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

    fun <- function(x) twoStepSplit(x, asNumeric = TRUE)
    saveMat(annotation[, c("PROBE", "CHRLOC")], fun = fun,
            pkgName = pkgName, pkgPath = pkgPath, envName = "CHRLOC")
    cat("CHRLOC complete\n")
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
        base <- tempfile()
        write.table(getProbe2SGD(yGenoUrl=getSrcUrl("YG")), base, sep = "\t", row.names = FALSE,
                    col.names = FALSE, quote = FALSE)
    }
    chrLengths <- getChrLengths("yeast")
    writeOrganism(pkgName, pkgPath, "yeast")
    writeChrLength(pkgName, pkgPath, chrLengths)
    repList <- getRepList("YEAST", srcObjs)
    repList[["PKGNAME"]] <- pkgName
    writeDocs(base, pkgName, pkgPath, version, author,
                      repList, "PKGNAME")
    if (lazyLoad)
      makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
}







