ABPkgBuilder <- function(baseName, srcUrls,
                         baseMapType=c("gb", "ug", "ll", "image", "refseq", "gbNRef"),
                         otherSrc=NULL,
                         pkgName, pkgPath, organism, version,
                         author, fromWeb = TRUE, lazyLoad = TRUE) {
    ## Build annotation data packages.
    ##
    ## pkgName - the name of the data package to be built (e.g. hgu95a)
    ##
    ## organism - a character string for the name of the organism of
    ##            concern (now can only be "human", "mouse", or "rat")
    ## pkgPath - a character string for the full path of an existing
    ##           directory where the built backage will be stored
    ## version - a character string for the version number
    ## baseName - a character string for the name of a file to be used as a
    ##            base file to base source data.
    ## baseMapType - a character string that can be either "gb", "ug",
    ##             "ll", "image", "refseq", or "gbNRef"  to
    ##             indicate whether the probe ids in baseName are mapped to
    ##             GenBack accession numbers, UniGene ids, LocsuLink ids,
    ##             Image clone ids, RefSeq ids, or a mixture of GenBank
    ##             accession numbers and RefSeq ids.
    ## otherSrc - a vector of named character strings for the names of files that
    ##          contain mappings between probe ids of baseName and
    ##          LobusLink ids that will be used to obtain the unified
    ##          mappings between probe ids of baseName and LocusLink ids
    ##          based on all the sources. The strings should not contain
    ##          any number.
    ## makeXML - a boolean to indicate whether an XML version will also be
    ##         generated.
    ## srcUrls - a vector of names character strings for the urls where
    ##         source data files will be retained. Valid sources are LocusLink,
    ##         UniGene, Golden Path, Gene Ontology, and KEGG. The names for the
    ##         character strings should be LL, UG, GP, GO, and KEGG,
    ##         respectively. LL and UG are required.
    ##
    ## author - a list with named elements "authors" containing a character
    ##          vector of author names and "maintainer" containing the complete
    ##          character string for the maintainer field, for example, "Jane
    ##          Doe <jdoe@doe.com>".
    ##
    ## This function remains to be long to avoid passing large data files
    ## around considering the size of the data files being processed.
    ##
    ## Copyright 2003, J. Zhang, all rights reserved.
    ##
    require("GO", quietly=TRUE) || stop("GO is needed to build data package")

    baseMapType <- match.arg(baseMapType)
    if(missing(srcUrls)){
         srcUrls <- getSrcUrl("all", organism = organism)
     }
     makeSrcInfo()

     srcObjs <- getSrcObjs(srcUrls, baseName, organism, baseMapType)
     if(!is.null(otherSrc) && !is.na(otherSrc))
     {
     	unified <- getUniMappings(baseName, srcObjs[["eg"]],
                               srcObjs[["ug"]], otherSrc, baseMapType)
     }
     else
     {
       unified <- getUniMappings(baseName, srcObjs[["eg"]],
                               srcObjs[["ug"]], NULL, baseMapType)
     }
     # Using the unified mapping as the base file, the data file
     # (ll_tmpl) can be parsed to get annoation data from LocusLink
     # using the correct parser
     #if(baseMapType == "ll"){
     #    parser(srcObjs[["ll"]]) <- file.path(.path.package("AnnBuilder"),
     #                            "scripts", "llParser4LL")
     #}else{
     #    parser(srcObjs[["ll"]]) <- file.path(.path.package("AnnBuilder"),
     #                            "scripts", "llParser")
     #}
     # Modify the setting for eg to be able to reuse it to parse other
     # Entrez Gene data sets
     ## Modify the setting for eg to be able to reuse it to parse other
     ## Entrez Gene data sets
     baseFile(srcObjs[["eg"]]) <- unified
     srcUrl(srcObjs[["eg"]]) <- srcUrls["EG"]
     fromWeb(srcObjs[["eg"]]) <- TRUE

     annotation <- getAnnData(srcObjs)
     
     createEmptyDPkg(pkgName = pkgName, pkgPath, force = TRUE)
     writeAnnData2Pkg(annotation, pkgName, pkgPath)
     revNames <- intersect(colnames(annotation), c("PMID", "PATH", "ENZYME"))
     if(length(revNames) != 0){
         writeReverseMap(annotation[, c("PROBE", revNames)],
                         pkgName, pkgPath)
     }
     if(!all(is.na(annotation[, "GO"]))){
         goCategories <- unlist(eapply(GOTERM, function(x) x@Ontology))
         goList <- getList4GO(goCategories,
                              sapply(annotation[, "GO"], twoStepSplit))
         names(goList) <- annotation[, "PROBE"]
         saveList(goList, pkgName, pkgPath, "GO")
         go2Probe <- mapGO2Probe(srcObjs[["eg"]], baseMapType)
         go2ProbeEnv <- saveMat(go2Probe, fun=twoStepSplit, pkgName=pkgName,
                                pkgPath=pkgPath, envName="GO2PROBE")
         go2AllEnv <- createGo2AllProbesEnv(go2ProbeEnv)
         lockEnvironment(go2AllEnv, bindings=TRUE)
         envName <- paste(pkgName, "GO2ALLPROBES", sep="")
         assign(envName, go2AllEnv)
         fName <- file.path(pkgPath, pkgName, "data", paste(envName, ".rda", sep=""))
         save(list=envName, file=fName, compress=TRUE)

         ## Create an empty SUMFUNC environment
         sumfunc <- as.list(rep(NA, length(annotation[, "PROBE"])))
         names(sumfunc) <- annotation[, "PROBE"]
         sumFuncEnv <- new.env(hash=TRUE, parent=NULL)
         sumFuncEnv <- l2e(sumfunc, sumFuncEnv)
         lockEnvironment(sumFuncEnv, bindings=TRUE)
         envName <- paste(pkgName, "SUMFUNC", sep="")
         assign(envName, sumFuncEnv)
         fName <- file.path(pkgPath, pkgName, "data",
                            paste(envName, ".rda", sep=""))
         save(list=envName, file=fName, compress=TRUE)
         
     }
     
     repList <- getRepList("all", resumeSrcUrl(srcObjs, organism))
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
     if (lazyLoad) {
       makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
     }
     ## clean up
     cleanSrcObjs(srcObjs)
}

resumeSrcUrl <- function(srcObjs, organism){
    if(!is.null(srcObjs[["ug"]])){
        srcUrl(srcObjs[["ug"]]) <- getSrcUrl("ug", organism = organism)
    }
    return(srcObjs)
}

# what should be "GO" - for GO package
#                "ALL" - for platform specific packages
#srcObjs a list of pubRepo objects
#
getRepList <- function(what, srcObjs){
    switch(toupper(what),
           GO = return(list(LLSOURCE = getRepSourceNBuilt("Entrez Gene:",
                                 srcObjs[["eg"]]),
                            GOSOURCE = getRepSourceNBuilt("Gene Ontology:",
                                 srcObjs[["go"]]), DATE = date())),
           ALL = return(list(LLSOURCE = getRepSourceNBuilt("Entrez Gene:",
                                        srcObjs[["eg"]]),
                             GOSOURCE = getRepSourceNBuilt("Gene Ontology:",
                                  srcObjs[["go"]]),
                             GPSOURCE = getRepSourceNBuilt("Golden Path:",
                                 srcObjs[["gp"]]),
                             KEGGSOURCE = getRepSourceNBuilt("KEGG:",
                                 srcObjs[["kegg"]]), DATE = date())),
           CHRLOC = return(list(GPSOURCE = getRepSourceNBuilt("Golden Path:",
                                 srcObjs[["gp"]]), DATE = date())),
           YEAST = return(list(LLSOURCE = getRepSourceNBuilt("Yeast Genome:",
                                   srcObjs[["yg"]]),
                               YGSOURCE = getRepSourceNBuilt("Yeast Genome:",
                                   srcObjs[["yg"]]),
                               GPSOURCE = getRepSourceNBuilt("Yeast Genome:",
                                   srcObjs[["yg"]]),
                               GOSOURCE = getRepSourceNBuilt("Gene Ontology:",
                                   srcObjs[["go"]]),
                               KEGGSOURCE = getRepSourceNBuilt("KEGG:",
                                   srcObjs[["kegg"]]), DATE = date())),
           LLMAP = return(list(LLSOURCE = getRepSourceNBuilt("Entrez Gene:",
                                 srcObjs[["eg"]]),
                          UGSOURCE = getRepSourceNBuilt("UniGene",
                          srcObjs[["ug"]]), DATE = date())),
           HG = return(list(HGSOURCE = getRepSourceNBuilt("HomoloGene:",
                                 srcObjs[["hg"]]), DATE = date())),
           KEGG = return(list( KEGGSOURCE = getRepSourceNBuilt("KEGG:",
                                   srcObjs[["kegg"]]), DATE = date())),
           UG = return(list(LLSOURCE = getRepSourceNBuilt("UniGene:",
                             srcObjs[["ug"]]), DATE = date())),
           stop(paste("Unknown argument", what)))
}

getRepSourceNBuilt <- function(name, object){
    if(!is.null(object)){
        return(paste(name, "\\\\url{", srcUrl(object), "}. Built: ",
                           escapeLatexChr(builtInfo(object)), sep = ""))
    }else{
        return(paste(name, "Not available"))
    }
}


# Map GO ids to probe ids that are directly associated with the GO ids
mapGO2Probe <- function(eg, baseMapType){
     if(baseMapType == "ll"){
            parser(eg) <- file.path(.path.package("AnnBuilder"),
                                    "scripts", "GO2ProbeParser4LL")
        }else{
            parser(eg) <- file.path(.path.package("AnnBuilder"),
                                    "scripts", "GO2ProbeParser")
        }
        options(show.error.messages = FALSE)
        go2Probe <- try(parseData(eg, eg@go, ncol = 3))
        options(show.error.messages = TRUE)
        if(inherits(go2Probe, "try-error")){
            stop(paste("Failed to parse LocusLink data because of:\n\n",
                       go2Probe))
        }
        if(is.null(go2Probe)){
            warning("No mappings between GO and probe ids found")
        }else if(nrow(go2Probe) == 1){
            go2Probe <- matrix(go2Probe[1:2], nrow = 1)
            colnames(go2Probe) <- c("GO", "GO2PROBE")
        }else{
            colnames(go2Probe) <- c("GO", "GO2PROBE", "COUNTS")
            # Drop column counts for now. Will get ride of the extra column
            go2Probe <- go2Probe[, c("GO", "GO2PROBE")]
            # Remove "" for GO id resulting from mapping probe ids to NA GO id
            go2Probe <- go2Probe[go2Probe[,1] != "", ]
        }
     return(go2Probe)
}


# Write reverse mappings between all the columns with column number >
# 1 to the first column of annData
writeReverseMap <- function(annData, pkgName, pkgPath){
    colNames <- colnames(annData)
    for(i in 2:length(colNames)){
        temp<- getReverseMapping(annData[, c(1, i)])
        if(length(temp) != 0){
            saveMat(temp, fun = splitEntry, pkgName = pkgName,
                    pkgPath = pkgPath, envName = paste(colNames[i],
                                       "2", colNames[1], sep = ""))
        }
    }
}

# Write mappings between the first column and each of the rest of the
# columns of annData
writeAnnData2Pkg <- function(annData, pkgName, pkgPath){
    #for (i in getUniColNames()) {
    #    env <- new.env(hash = TRUE, parent = NULL)
    #    if (i == "LOCUSID") {
    #        multiassign(annData[, "PROBE"], as.integer(annData[,
    #            i]), env)
    #    }
    #    else {
    #        multiassign(annotation[, "PROBE"], as.vector(annData[,
    #            i]), env)
    #    }
    #    assign(paste(pkgName, i, sep = ""), env)
    #    save(list = paste(pkgName, i, sep = ""), file = file.path(pkgPath,
    #        pkgName, "data", paste(pkgName, i, ".rda", sep = "")))
    #}
    #for (i in intersect(colnames(annData), getMultiColNames())) {
    #    env <- new.env(hash = TRUE, parent = NULL)
    #    multiassign(annData[, "PROBE"], lapply(annData[,
    #        i], splitEntry), env)
    #    assign(paste(pkgName, i, sep = ""), env)
    #    save(list = paste(pkgName, i, sep = ""), file = file.path(pkgPath,
    #        pkgName, "data", paste(pkgName, i, ".rda", sep = "")))
    #}
    #for (i in intersect(colnames(annData), c("CHRLOC"))) {
    #    env <- new.env(hash = TRUE, parent = NULL)
    #    if (i == "CHRLOC") {
    #        multiassign(annData[, "PROBE"], lapply(annData[,
    #            i], twoStepSplit, asNumeric = TRUE), env)
    #    }
    #    else {
    #        multiassign(annData[, "PROBE"], lapply(annData[,
    #            i], twoStepSplit), env)
    #    }
    #    assign(paste(pkgName, i, sep = ""), env)
    #    save(list = paste(pkgName, i, sep = ""), file = file.path(pkgPath,
    #        pkgName, "data", paste(pkgName, i, ".rda", sep = "")))
    #}

    colNames <- colnames(annData)
    colNames <- colNames[!colNames %in% c("PROBE", "GO")]

    # Write data to the package for one to one mappings
    for(i in colNames){
        switch(i,
               LOCUSID = fun <- function(x) as.integer(x),
               CHRLOC = fun <- function(x) twoStepSplit(x,
                                              asNumeric = TRUE),
               fun <- function(x) splitEntry(x))
        saveMat(annData[, c("PROBE", i)], fun = fun,
                pkgName = pkgName, pkgPath = pkgPath, envName = i)
    }
}

getAnnData <- function(srcObjs){
    annotation <- read.delim(baseFile(srcObjs[["eg"]]), sep = "\t",
                             header = FALSE, as.is = TRUE)
    annotation <- annotation[, 1:3]
    colnames(annotation) <- c("PROBE", "ACCNUM", "LOCUSID")
    options(show.error.messages = FALSE)
    # Parse gene2go.gz
    parser(srcObjs[["eg"]]) <- getBaseParsers("eggo") 
    go <- try(parseData(srcObjs[["eg"]], srcObjs[["eg"]]@go, 
                                ncol = 2, mergeKey = FALSE))
    colnames(go) <- c("PROBE", "GO")
    options(show.error.messages = TRUE)
    if(inherits(annotation, "try-error")){
        stop(paste("Parsing Entrez Gene gene2go.gz failed because of:\n\n",
                   annotation))
    }
    if(nrow(go) > 0)
      annotation <- merge(annotation, go, by = "PROBE", all.x = TRUE)
    options(show.error.messages = FALSE)
    # Parse gene2pubmid.gz
    parser(srcObjs[["eg"]]) <- getBaseParsers("egpubmed") 
    pubmed <- try(parseData(srcObjs[["eg"]], srcObjs[["eg"]]@pubmed, 
                                ncol = 2, mergeKey = FALSE))
    colnames(pubmed) <- c("PROBE", "PMID")
    options(show.error.messages = TRUE)
    if(inherits(annotation, "try-error")){
        stop(paste("Parsing Entrez Gene gene2pubmed.gz failed because of:\n\n",
                   annotation))
    }
    if(nrow(pubmed) > 0)
      annotation <- merge(annotation, pubmed, by = "PROBE", all.x = TRUE)
    # Parse gene_info.gz
    parser(srcObjs[["eg"]]) <- getBaseParsers("eginfo") 
    geneinfo <- try(parseData(srcObjs[["eg"]], srcObjs[["eg"]]@info, 
                                ncol = 5, mergeKey = FALSE))
    colnames(geneinfo) <- c("PROBE", "SYMBOL", "GENENAME", "CHR", "MAP")
    options(show.error.messages = TRUE)
    if(inherits(annotation, "try-error")){
        stop(paste("Parsing Entrez Gene gene_info.gz failed because of:\n\n",
                   annotation))
    }
    if(nrow(geneinfo) > 0)
      annotation <- merge(annotation, geneinfo, by = "PROBE", all.x = TRUE)
    # Parse gene2refseq.gz
    parser(srcObjs[["eg"]]) <- getBaseParsers("egrefseq") 
    refseq <- try(parseData(srcObjs[["eg"]], srcObjs[["eg"]]@refseq, 
                                ncol = 2, mergeKey = FALSE))
    colnames(refseq) <- c("PROBE", "REFSEQ")
    options(show.error.messages = TRUE)
    if(inherits(annotation, "try-error")){
        stop(paste("Parsing Entrez Gene gene2refseq.gz failed because of:\n\n",
                   annotation))
    }
    if(nrow(refseq) > 0)
      annotation <- merge(annotation, refseq, by = "PROBE", all.x = TRUE)
    # Parse gene2unigene.gz
    options(show.error.messages = FALSE)
    parser(srcObjs[["eg"]]) <- getBaseParsers("egunigene") 
    unigene <- try(parseData(srcObjs[["eg"]], srcObjs[["eg"]]@unigene, 
                                ncol = 2, mergeKey = FALSE))
    colnames(unigene) <- c("PROBE", "UNIGENE")
    options(show.error.messages = TRUE)
    if(inherits(annotation, "try-error")){
        stop(paste("Parsing Entrez Gene gene2unigene failed because of:\n\n",
                   annotation))
    }
    if(nrow(unigene) > 0)
      annotation <- merge(annotation, unigene, by = "PROBE", all.x = TRUE)
    # Parse mm2gene.gz
    parser(srcObjs[["eg"]]) <- getBaseParsers("egmim") 
    mim <- try(parseData(srcObjs[["eg"]], srcObjs[["eg"]]@mim, 
                                ncol = 2, mergeKey = FALSE))
    colnames(mim) <- c("PROBE", "OMIM")
    options(show.error.messages = TRUE)
    if(inherits(annotation, "try-error")){
        stop(paste("Parsing Entrez Gene mim2gene failed because of:\n\n",
                   annotation))
    }
    if(nrow(mim) > 0)
      annotation <- merge(annotation, mim, by = "PROBE", all.x = TRUE)

    if(!is.null(srcObjs[["gp"]])){
    # Get strand data by processing data from GoldenPath
        options(show.error.messages = FALSE)
        strand <- try(getStrand(srcObjs[["gp"]]))
        options(show.error.messages = TRUE)
        if(inherits(strand, "try-error")){
            warning(paste("Failed to parse Golden Path data because of:\n\n",
                       strand))
        }else{
            annotation <- merge(annotation, strand, by = "LOCUSID",
                                all.x = TRUE)
        }
    }
    if(!is.null(srcObjs[["kegg"]])){
        # Adding KEGG pathway data to annoataion
        options(show.error.messages = FALSE)
        pathNEnzyme <- try(mapLL2ECNPName(srcObjs[["kegg"]]))
        options(show.error.messages = TRUE)
        if(inherits(pathNEnzyme, "try-error")){
            warning("Failed to include KEGG data")
        }else{
            print(paste(dim(pathNEnzyme$llpathname)[[1]], dim(pathNEnzyme$llpathname)[[2]], length(pathNEnzyme), sep=" "))
            if(dim(pathNEnzyme$llpathname)[[1]] > 0)
            {
               annotation <- merge(annotation, pathNEnzyme$llpathname,
                                by = "LOCUSID", all.x = TRUE)
               annotation <- merge(annotation, pathNEnzyme$llec,
                                by = "LOCUSID", all.x = TRUE)
            }
            else
            {
               PATH <- rep(NA, length(annotation[, "PROBE"]))
	       ENZYME <- rep(NA, length(annotation[, "PROBE"]))
               GO <- rep(NA, length(annotation[, "PROBE"]))
               annotation <- data.frame(annotation, PATH, ENZYME, GO)
            }
        }
    }
    annotation <- as.matrix(annotation)
    # Convert "" to NA
    annotation[annotation == ""] <- NA
    # Probe id can not be NULL
    annotation <- annotation[!is.null(annotation[, "PROBE"]),]

    return(annotation)
}

.parseEGUniGene <- function(eg){
  options(show.error.messages = FALSE)
  unigene <- try(parseData(eg, eg@unigene, ncol = 2, mergeKey = FALSE))
  colnames(unigene) <- c("PROBE", "UNIGENE")
  options(show.error.messages = TRUE)
  if(inherits(unigene, "try-error")){
    stop(paste("Parsing Entrez Gene gene2unigene failed because of:\n\n",
               unigene))
  }
  return(unigene)
}
  
unifyMappings <- function(base, eg, ug, otherSrc){

    trusted <- NULL
    # Get the unified mappings betwee probe ids and locusLink ids
    # based on multiple sources.
    if(!is.null(eg)){
        options(show.error.messages = FALSE)
        llMapping <- try(parseData(eg, eg@accession))
        options(show.error.messages = TRUE)
        if(inherits(llMapping, "try-error")){
            warning(paste("Failed to get or parse Entrez Gene data",
                       "because of:\n\n",llMapping))
        }else{
            colnames(llMapping) <- c("PROBE", "LOCUSLINK")
        }
    }
    if(!is.null(ug) && !is.na(srcUrl(ug))){
        options(show.error.messages = FALSE)
        ugMapping <- try(parseData(ug))
        options(show.error.messages = TRUE)
        if(inherits(ugMapping, "try-error")){
            warning(paste("Failed to get or parse UniGene data becaus",
                          "of:\n\n", ugMapping))
        }else{
            colnames(ugMapping) <- c("PROBE", "UNIGENE")
        }
    }
    # Merge data from all the sources by probe id based on the base
    # file. The probe id is the first column ("V1") in this case.
    # Merge two at a time.
    if(!is.null(eg) && nrow(llMapping) > 0){
        merged <- merge(base, llMapping, by = "PROBE", all.x = TRUE)
        trusted <- "LOCUSLINK"
         # Free the space
        llMapping <- NULL
    }else{
        merged <- base
    }
    if(!is.null(ug) && nrow(ugMapping) > 0){
        merged <- merge(merged, ugMapping, by = "PROBE", all.x = TRUE)
        trusted <- c(trusted, "UNIGENE")
        ugMapping <- NULL
    }
    # The following three mappings are for Affymetrix gene chips.
    # Skip if the target of annotation is not a Affymetrix chip with
    # existing mapping data file in the data directory of AnnBuilder.
    if(!is.null(otherSrc) && length(otherSrc) > 0){
        names(otherSrc) <- paste("OTHER", names(otherSrc))
        for(i in names(otherSrc)){
            options(show.error.messages = FALSE)
            temp <- try(matrix(scan(otherSrc[i], what = "character",
                         sep ="\t", quote = "", quiet = TRUE,
                         strip.white = TRUE), ncol = 2, byrow = TRUE))
            options(show.error.messages = TRUE)
            if(inherits(temp, "try-error")){
                stop(paste("File", otherSrc[i],
                           "may not be valid or have two columns"))
            }
            colnames(temp) <- c("PROBE", i)
            merged <- merge(merged, temp, by = "PROBE", all.x = TRUE)
        }
    }
    # Change "NA" tp NA
    merged <- as.matrix(merged)
    merged[merged == "NA"] <- NA
    # Get the unified mapping based on all the sources used
    unified <- resolveMaps(merged, trusted = trusted,
                            srcs = c(names(otherSrc)))
    return(unified)
}

getUniMappings <- function(baseName, eg, ug, otherSrc, baseMapType){

    base <- getBaseFile(baseName)

    if(baseMapType == "refseq"){
        unified <- unifyMappings(base, eg, NULL, otherSrc)
    }else if(baseMapType != "ll"){
        #ug <- UG(srcUrl = srcUrl(ug),
        #         parser = parser(ug), baseFile = baseName(ug),
        #         organism = organism)
        if(baseMapType == "image"){
            unified <- unifyMappings(base, NULL, ug, otherSrc)
        }else{
            unified <- unifyMappings(base, eg, ug, otherSrc)
        }
    }else{
        unified <- baseName
    }
    return(unified)
}

#mapImage2Ref <- function(baseName, ll, ug, organism){
#    base <- getBaseFile(baseName)
#    f <- loadFromUrl(paste(getSrcUrl("gp", organism = organism),
#                      "imageClone.txt.gz", sep = ""), tempdir())
#    data <- matrix(scan(f, what = "", sep = "\t", strip.white = TRUE,
#                        quiet = TRUE), byrow = TRUE, ncol = 4)[, 1:2]
#    colnames(data) <- c("ACC", "REFNGB")

#    merged <-mergeRowByKey(merge(base, data, by = "ACC", all.x = TRUE))
    #merged <- mergeRowByKey(merged)


#}

getBaseFile <- function(baseName){

    path <- file.path(.path.package("AnnBuilder"), "data")
    # Read in the base file
    options(show.error.messages = FALSE)
    base <- try(matrix(scan(baseName, what = "", sep = "\t", quote = "",
                        quiet = TRUE, strip.white = TRUE), ncol = 2,
                       byrow = TRUE))
    options(show.error.messages = TRUE)
    if(inherits(base, "try-error")){
        stop(paste("Base file", baseName,
                   "may not be valid or have two columns"))
    }
    colnames(base) <- c("PROBE", "ACC")

    return(base)
}


getSrcObjs <- function(srcUrls, baseName, organism,
                  baseMapType = c("gb", "ug", "ll", "image", "refseq",
                  "gbNRef"), fromWeb = TRUE){

    baseMapType <- match.arg(baseMapType)
    switch(baseMapType,
           gb = baseParser <- getBaseParsers("gb"),
           ug = baseParser <- getBaseParsers("ug"),
           image = baseParser <- getBaseParsers("image"),
           ll = baseParser <- getBaseParsers("ll"),
           refseq = baseParser <- getBaseParsers("refseq"),
           gbNRef = baseParser <- getBaseParsers("gbNRef"),
           stop("Invalid imput for baseMapType"))

    srcObjs <- list()
    for(i in names(srcUrls)){
       switch(toupper(i),
              EG = srcObjs[["eg"]] <- EG(loadFromUrl(paste(srcUrls[i],
                              getEGAccName(), sep = "/")),
                              parser = baseParser["EG"], baseFile = baseName,
                              fromWeb = FALSE),
              UG =  srcObjs[["ug"]] <- UG(loadFromUrl(srcUrls[i]),
                              parser = baseParser["UG"], baseFile = baseName,
                              organism = organism, fromWeb = FALSE),
              KEGG = srcObjs[["kegg"]] <- KEGG(srcUrl = srcUrls[i],
                            organism = organism, built = getSrcBuilt("kegg"),
                            fromWeb = TRUE),
              GO =  srcObjs[["go"]] <- GO(srcUrl = srcUrls[i],
                          built = getSrcBuilt("go"), fromWeb = TRUE),
              GP = srcObjs[["gp"]] <- GP(srcUrl = srcUrls[i],
                          organism = organism,
                          built = getSrcBuilt("gp", organism = organism),
                          fromWeb = TRUE))
    }

    return(srcObjs)
}

getEGAccName <- function(){
  return("gene2accession.gz")
}


# If the source file is not from a web site, the file has been stored
# locally and should be droped when done using it
cleanSrcObjs <- function(srcObjs){
    for(i in names(srcObjs)){
        if(!fromWeb(srcObjs[[i]])){
            unlink(srcUrl(srcObjs[[i]]))
        }
    }
}

getBaseParsers <- function(baseMapType = c("gb", "ug", "image", "ll",
    "refseq", "gbNRef", "ll2gb", "gb2ll", "eggo", "eginfo", "egrefseq",
                             "egpubmed", "egunigene", "egmim")){
    type <- match.arg(baseMapType)
    path <- file.path(.path.package("AnnBuilder"), "scripts")
    switch(type,
           gb = return(c(EG = file.path(path, "egAccParser"),
                         UG = file.path(path, "gbUGParser"))),
           ug = return(c(EG = file.path(path, "egUGParser"),
                         UG = file.path(path, "ugUGParser"))),
           image = return(c(UG = file.path(path, "imageUGParser"),
                          LL = "NA")),
           ll = return(c(LL = file.path(path, "llParser"),
                         UG = "NA")),
           refseq = return(c(EG = file.path(path, "egAccParser"),
                             UG = "NA")),
           gbNRef = return(c(EG = file.path(path, "egAccParser"),
                             UG = file.path(path, "gbUGParser"))),
           ll2gb = return(c(UG = file.path(path, "ll2GBUGParser"))),
           gb2ll = return(c(UG = file.path(path, "gb2LLUGParser"))),
           eggo = return(EG = file.path(path, "egGoParser")),
           eginfo = return(EG = file.path(path, "egInfoParser")),
           egrefseq = return(EG = file.path(path, "egRefseqParser")),
           egpubmed = return(EG = file.path(path, "egPubMedParser")),
           egunigene = return(EG = file.path(path, "egUnigeneParser")),
           egmim = return(EG = file.path(path, "egMimParser")),
           stop("Invalid type"))
}

createEmptyDPkg <- function(pkgName, pkgPath,
                            folders = c("man", "R", "data"), force = TRUE){
    if(file.exists(file.path(pkgPath, pkgName))){
        if(!force){
            stop(paste("Package", pkgName, "already exists"))
        }else{
            unlink(file.path(pkgPath, pkgName), TRUE)
        }
    }
    dir.create(file.path(pkgPath, pkgName))
    for(i in folders){
        dir.create(file.path(pkgPath, pkgName, i))
    }
}
# Split multiple entry for a given mapping
splitEntry <- function(dataRow, sep = ";", asNumeric = FALSE){
    if(is.na(dataRow) || is.null(dataRow) || dataRow == ""){
        return(NA)
    }else{
        if(asNumeric){
            return(unique(as.numeric(unlist(strsplit(dataRow, sep)))))
        }else{
            return(unique(unlist(strsplit(dataRow, sep))))
        }
    }
}
# Split multiple entry with two separaters (e. g. 12345@18;67891@18)
twoStepSplit <- function(dataRow, entrySep = ";", eleSep = "@",
    asNumeric = FALSE){
    splitEle <- function(entry){
        if(!is.na(entry)){
            temp <- unlist(strsplit(entry, eleSep), use.names = FALSE)
            if(asNumeric){
                elements <- as.numeric(temp[1])
            }else{
                elements <- temp[1]
            }
            names(elements) <- temp[2]
            return(elements)
        }else{
            return(NA)
        }
    }
    temp <- unique(unlist(strsplit(dataRow, entrySep), use.names = FALSE))
    # Force "NA" to be NA
    temp[temp == "NA"] <- NA
    return(unlist(lapply(temp, splitEle), use.names = TRUE))
}

# This function takes two data columns from a data file (e. g. two
# columns) and generates an environment object with
# values in one of the columns as keys and values in the other as
# values. A given row in either column may have multiple values
# separated by a separator(e. g. "a;b") and has to be dealt with.
#
# Copyright 2002, J. Zhang, all rights reserved
#

cols2Env <- function(cols, colNames, keyColName = colNames[1], sep = ";"){

    # Environment to return
    dataEnv <- new.env(hash = TRUE, parent = NULL)

    # Force a matrix
    cols <- as.matrix(cols)
    if(!missing(colNames)){
        colnames(cols) <- colNames
    }

    if(ncol(cols) != 2){
        stop("The input data do not have two columns")
    }
    # Determine the column name
    valColName <- colNames[colNames != keyColName]
    # Match cols by finding matches for all the potential multiple
    # matches
    cols <- matchAll(cols, keyColName)
    colnames(cols) <- colNames
    #Get unduplicated row
    unDups <- cols[!duplicated(cols[,keyColName]),]
    colnames(unDups) <- colNames
    # Write unique keys and values
    multiassign(unDups[, keyColName], unDups[, valColName], dataEnv)
    # Process the duplications
    dups <- cols[duplicated(cols[,keyColName]),]
    colnames(dups) <- colNames
    if(nrow(dups) > 0){
        for(i in 1:nrow(dups)){
            assign(dups[i,keyColName], unique(c(get(dups[i, keyColName],
                                   dataEnv), dups[i,valColName])), dataEnv)
        }
    }

    return(dataEnv)
}

matchAll <- function(cols, keyColName){
    matched <- NULL

    for(i in 1:nrow(cols)){
        temp <- matchOneRow(cols[i,], keyColName)
        if(!is.null(temp)){
            if(is.null(matched)){
                matched <- temp
            }else{
                matched <- rbind(matched, temp)
            }
        }
    }
    return(matched)
}

matchOneRow <- function(cols, keyColName, sep = ";"){

    doMatch <- function(col1id){
        if(is.null(matched)){
            matched <<- cbind(col1id, col2)
        }else{
            matched <<- rbind(matched, cbind(col1id, col2))
        }
    }

    matched <- NULL
    # Key can not be any of these
    if(cols[keyColName] == "" || is.na(cols[keyColName])||
       is.null(cols[keyColName])){
        return(NULL)
    }

    col1 <- unlist(strsplit(cols[1], sep), use.names = FALSE)
    col2 <- unlist(strsplit(cols[2], sep), use.names = FALSE)
    # If no mapping is found for key, assign NA to key
    if(match(keyColName, names(cols)) == 1 && length(col2) == 0){
        return(cbind(col1, NA))
    }else if(match(keyColName, names(cols)) == 2 && length(col1) == 0){
        return(cbind(NA, col2))
    }
    # Otherwise, do a full matching
    sapply(col1, doMatch)

    return(matched)
}


getDirContent <- function(dirName, exclude = NULL){
    if(is.null(exclude)){
        return(list.files(dirName, full.names = TRUE))
    }else{
        temp <- list.files(dirName, full.names = TRUE)
        index <- unlist(lapply(exclude, grep, temp))
        if(length(index) > 0){
            return(temp[-index])
        }else{
            return(temp)
        }
    }
}

# Returns a vector for the names of columns that have multiple entries
# separated by a ";" in some of the cells
getMultiColNames <- function(){
     return(c("PMID", "PATH", "ENZYME", "CHR", "UNIGENE", "HGID",
              "OMIM", "REFSEQ"))
}

getUniColNames <- function(){
    return(c("ACCNUM", "LOCUSID", "GENENAME", "SYMBOL", "MAP",
             "GRIF", "SUMFUNC"))
}

getTypeColNames <- function(){
    return(c("GO", "CHRLOC"))
}

# goNCat = a named vector with GO category as the values and GO id as
#          the names
# goNEvi = a list of named vectors with GO ids as values for vectors
#          and evidence code as names for vector values
getList4GO <- function(goNCat, goNEvi){
    procOne <- function(goids){
        if(is.null(goids) || is.na(goids)){
            return(NA)
        }else{
            temp <- cbind(goids, names(goids),goNCat[goids])
            rownames(temp) <- goids
            return(apply(temp, 1, vect2List,
                            vectNames = c("GOID", "Evidence", "Ontology")))
        }
    }
    temp <- sapply(goNEvi, procOne)
    # Names from sapply can be very long
    names(temp) <- 1:length(temp)
    return(temp)
}


vect2List <- function(vector, vectNames){
    tempList <- as.list(vector)
    names(tempList) <- vectNames
    return(tempList)
}


# This function writes the organism data
writeOrganism <- function(pkgName, pkgPath, organism){
    assign(paste(pkgName, "ORGANISM", sep = ""), organism)
    save(list = paste(pkgName, "ORGANISM", sep = ""),
             file =  file.path(pkgPath, pkgName, "data",
             paste(pkgName, "ORGANISM.rda", sep = "")))
}

# This function writes chromosome length data
writeChrLength <- function(pkgName, pkgPath, chrLengths){
    assign(paste(pkgName, "CHRLENGTHS", sep = ""), chrLengths)
    save(list = paste(pkgName, "CHRLENGTHS", sep = ""),
             file =  file.path(pkgPath, pkgName, "data",
             paste(pkgName, "CHRLENGTHS.rda", sep = "")))
}

# This function estimates the total lengths of chromosomes by finding
# out the maximum chromosome loaction and then increasing the number by 1000
findChrLength <- function(organism, srcUrl = getSrcUrl("GP", organism)){
    chrLengths <- vector()

    estimateLength <- function(chroNum){
        locs <- locatData[locatData[,1] == as.character(chroNum),]
        if(nrow(locs) == 0){
            temp <- NA
        }else if (nrow(locs) == 1){
            temp <- as.numberic(locs[, 2])
        }else {
            temp <- max(as.numeric(locs[, 2]))
        }
        chrLengths[as.character(chroNum)] <<- temp
    }
    locatData <- getGPData(paste(srcUrl, "refGene.txt.gz", sep = ""),
                           ncol = 10, keep = c(2,4))
    # Remove "chr" from chromosome numbers
    locatData[,1] <- gsub("^chr", "\\", locatData[, 1])
    chromosomes <- unique(locatData[,1])
    # Do numbered chromosomes first
    options(warn = -1)
    numbered <- chromosomes[!is.na(as.numeric(chromosomes))]
    options(warn = 0)
    tt <- sapply(sort(as.numeric(numbered)), estimateLength)
    # Do others
    switch(toupper(organism),
           "MUS MUSCULUS" = ,
           "RATTUS NORVEGICUS" = ,
           "HOMO SAPIENS" = sapply(c("X", "Y"), estimateLength),
           stop(paste("Do not know what to do with", organism)))
    # Remove entries with random locations
    # Add 1000 bases
    return(chrLengths + 1000)
}


getChrLengths <- function(organism){
    switch(toupper(organism),
           "HOMO SAPIENS" = return(getHumanChrLengths()),
           "MUS MUSCULUS" = return(getMouseChrLengths()),
           "RATTUS NORVEGICUS" = return(getRatChrLengths()),
           YEAST = return(getYeastChrLengths()),
           stop("Unknown organism"))
}

getHumanChrLengths <- function(){
    return(c("1" = 246127941, "2" = 243615958,
             "3" = 199344050, "4" = 191731959, "5" = 181034922,
             "6" = 170914576, "7" = 158545518, "8" = 146308819,
             "9" = 136372045, "10" = 135037215, "11" = 134482954,
             "12" = 132078379, "13" = 113042980, "14" = 105311216,
             "15" = 100256656, "16" = 90041932, "17" = 81860266,
             "18" = 76115139, "19" = 63811651, "20" = 63741868,
             "21" = 46976097, "22" = 49396972, "X" =153692391,
             "Y" = 50286555, "M" = 16571))
}

getMouseChrLengths <- function(){
    return(c("1" = 195869683, "2" = 181423755, "3" = 160674399,
             "4" = 152921959, "5" = 149719773, "6" = 149950539,
             "7" = 134401573, "8" = 128923138, "9" = 124467299,
             "10" = 130738012, "11" = 122862689, "12" = 114462600,
             "13" = 116242670, "14" = 115844145, "15" = 104111694,
             "16" = 98986639, "17" = 93529596, "18" = 91041441,
             "19" = 61093376, "X" = 149996094, "Y" = 687262))
}

getRatChrLengths <- function(){
    return(c("1" = 268121971, "2" = 258222147, "3" = 170969371,
             "4" = 187371129, "5" = 173106704, "6" = 147642806,
             "7" = 143082968, "8" = 129061546, "9" = 113649943,
             "10" = 110733352, "11" = 87800381, "12" = 46649226,
             "13" = 111348958, "14" = 112220682, "15" = 109774626,
             "16" = 90224819, "17" = 97307196, "18" = 87338544,
             "19" = 59223525, "20" = 55296979, "X" = 160775580,
             "Y" = NA))
}

getYeastChrLengths <- function(){
    return(c("1" = 230210, "2" = 813138, "3" =316613,
             "4" = 1531914, "5" = 576869, "6" = 270148,
             "7" = 1090944, "8" = 562639, "9" = 439885,
             "10" = 745446, "11" = 666445, "12" = 1078173,
             "13" = 924430, "14" =784328 , "15" =1091285,
             "16" = 948060, "MT" = 85779))
}

nameGOByCat <- function(GOWithEvi, goCat){
    if(is.na(GOWithEvi)){
        return(NA)
    }
    goids <- twoStepSplit(GOWithEvi)
    goids <- goCat[is.element(goCat[, 1], goids),]

    if(is.null(nrow(goids))){
        if(length(goids) > 0){
            return(paste(goids, sep = "", collapse = "@"))
        }else{
            return(GOWithEnv)
        }
    }else{
        return(paste(goids[,1], goids[,2], sep = "@", collapse = ";"))
    }
}

saveList <- function(dList, pkgName, pkgPath, envName){
    env <- new.env(hash = TRUE, parent = NULL)
    if(length(dList) != 0){
        l2e(dList, env)
        #multiassign(names(dList), dList, env)
    }
    lockEnvironment(env, bindings = TRUE)
    assign(paste(pkgName, envName, sep = ""), env)
    save(list = paste(pkgName, envName, sep = ""), file = file.path(pkgPath,
        pkgName, "data", paste(pkgName, envName, ".rda", sep = "")))
}

saveMat <- function(data, pkgName, pkgPath, envName, keyCol = 1,
                         valCol = 2, fun = function(x) x){
    env <- new.env(hash = TRUE, parent = NULL)
    if(is.null(nrow(data))){
        if(length(data) != 0){
            assign(data[keyCol], fun(data[valCol]), env)
        }
    }else{
        if(length(data) != 0){
            multiassign(data[,keyCol], lapply(as.vector(data[, valCol]),
                                              fun), env)
        }
    }
    lockEnvironment(env, bindings = TRUE)
    assign(paste(pkgName, envName, sep = ""), env)
    save(list = paste(pkgName, envName, sep = ""), file = file.path(pkgPath,
        pkgName, "data", paste(pkgName, envName, ".rda", sep = "")))
    env
}

# Create a datalist in the data directory so that lazy loading works
writeDatalist <- function(pkgName, pkgPath){

    dFiles <- gsub("\\..*", "",
                   list.files(file.path(pkgPath, pkgName, "data")))
    write.table(matrix(dFiles, ncol = 1),
                file = file.path(pkgPath, pkgName, "data", "datalist"),
                quote = FALSE, col.names = FALSE, row.names = FALSE)
}


addGoAllProbeData <- function(go2ProbeEnv, goOffspringEnv, allProbesEnv) {
    ## Add GO 2 all probe data to the allProbesEnv.  Specifially, this
    ## environment has keys of goIds and the value is a vector of probeIds such
    ## that the probeIds are associated with at least one of this goIds
    ## offspring.
    ##
    ## What we do is for each goId in the goOffspringEnv, we look up any matches
    ## in go2ProbeEnv and add these to allProbesEnv.
    ##
    for (goId in ls(goOffspringEnv)) {
        offspring <- c(goId, goOffspringEnv[[goId]])
        offspring <- offspring[!is.na(offspring)]
        ansList <- mget(offspring, go2ProbeEnv, ifnotfound=NA)
        ansList <- ansList[!is.na(ansList)]
        ans <- unlist(ansList, use.names=FALSE)
        names(ans) <- unlist(lapply(ansList, function(x) names(x)), use.names=FALSE)
        if (length(ans) == 0)
          next
        allProbesEnv[[goId]] <- ans
    }
    allProbesEnv
}
        

createGo2AllProbesEnv <- function(go2ProbeEnv) {
    ## assumes GO package is loaded
    allProbesEnv <- new.env(hash=TRUE, parent=NULL)
    addGoAllProbeData(go2ProbeEnv, GOCCOFFSPRING, allProbesEnv)
    addGoAllProbeData(go2ProbeEnv, GOMFOFFSPRING, allProbesEnv)
    addGoAllProbeData(go2ProbeEnv, GOBPOFFSPRING, allProbesEnv)
    allProbesEnv
}

