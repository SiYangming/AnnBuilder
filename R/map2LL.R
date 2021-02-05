# Makes the orgNameLLMappings packages

map2LL <- function(pkgName, pkgPath, organism, version, author,
                   eg = EG(parser = file.path(.path.package("AnnBuilder"),
                             "scripts", "egLLMappingUGParser")),
                   lazyload = TRUE){
    # Write AnnInfo to .GlobalEnv to make Annotation information available
    makeSrcInfo()
    # Create a data package with no data
    createEmptyDPkg(pkgName, pkgPath, force = TRUE)
    # Read data with mappings between Gene id and UG
    
    gene2ug = loadFromUrl(paste(getSrcUrl("eg"), eg@unigene, sep="/"))
    
    ug <- read.table(gene2ug, na.strings="-", sep="\t", quote="", colClasses=c("numeric","character"), comment.char="#", col.names=c("PROBE","UNIGENE"), header=FALSE)
    
    # Retain values for the organism
    ug <- ug[grep(getOrgName(organism, "scientific"), ug[,2]), ]
    
    
    ll2UG <- split(as.character(ug[,2]), as.character(ug[,1]))
    saveList(ll2UG, pkgName = pkgName, pkgPath = pkgPath, envName = "LL2UG")
    
    ug2LL <- split(as.numeric(ug[,1]), as.character(ug[,2]))   
    saveList(ug2LL, pkgName = pkgName, pkgPath = pkgPath, envName = "UG2LL")
    
    
    # Map LL 2 Acc number
    ugurl <- loadFromUrl(getSrcUrl("ug", organism))
    ugObj <- UG(srcUrl = ugurl, organism = organism, parser = "",
                                    fromWeb = FALSE)
    egurl <- loadFromUrl(paste(srcUrl(eg), eg@accession, sep = "/"))
    mapLLNGB(organism, pkgName, pkgPath, ugUrl = ugurl, egUrl = egurl,
                   fromWeb = FALSE)
    # Parse LL_tmpl data to get the mappings between LL and PMID
    oriParser <- file.path(.path.package("AnnBuilder"), "templates",
                            "llMappingPmid.temp")
    newParser <- tempfile("llMappingParser")

    # Write organism to the parser to be used
    orgCode <- getOrgNameNCode()
    copySubstitute(oriParser, newParser,
                   list(ORGCODE = names(orgCode[orgCode == organism])), "#")
    parser(eg) <- newParser
    options(show.error.messages = FALSE)
    pubmed <- try(parseData(eg, eg@pubmed, ncol = 2))
    options(show.error.messages = TRUE)
    if(inherits(pubmed, "try-error")){
        stop(paste("Failed to get or parse gene2pubmed.gz  becaus of:\n\n",
                   putmed))
    }
    # Do pmid first
    saveData2Env(pubmed, pkgName = pkgName,
                 pkgPath = pkgPath, envName = "LL2PMID")
    # Now do the reverse mapping for PMID
    pmid2LL <- getReverseMapping(pubmed)
    saveData2Env(pmid2LL, pkgName = pkgName,
                         pkgPath = pkgPath, envName = "PMID2LL")
    # Do GO ids
    oriParser <- file.path(.path.package("AnnBuilder"), "templates",
                            "llMappingGo.temp")
    copySubstitute(oriParser, newParser,
                   list(ORGCODE = names(orgCode[orgCode == organism])), "#")
    parser(eg) <- newParser
    options(show.error.messages = FALSE)
    go <- try(parseData(eg, eg@go, ncol = 2))
    options(show.error.messages = TRUE)
    if(inherits(go, "try-error")){
        stop(paste("Failed to get or parse gene2go.gz becaus of:\n\n",
                   go))
    }
    saveData2Env(go, pkgName = pkgName,
                         pkgPath = pkgPath, envName = "LL2GO")
    # Do the reverse mapping for go
    go2LL <- reverseMap4GO(go, sep = ";")
    saveData2Env(go2LL, twoStepSplit, pkgName, pkgPath, "GO2LL")
    # Write man pages and so on
    repList <- getRepList("LLMAP", list(ll = eg, ug = ugObj))
    repList[["PKGNAME"]] <- pkgName
    writeDocs("", pkgName, pkgPath, version, author, repList, "PKGNAME")
    writeDatalist(pkgName, pkgPath)
    if(lazyload){
           makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
     }
    unlink(c(ugurl, egurl, newParser))
    return(invisible())
}

#getUrl4Org <- function(organism){
#    switch(toupper(organism),
#           HUMAN = return(),
#           MOUSE = return(),
#           RAT = return(),
#           stop("Organism not supported"))
#}

getOrgName <- function(organism, what = c("common", "scientific")){
    what <- tolower(match.arg(what))
    names <- cbind(c("human", "mouse", "rat", "fly"),
                   c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Drosophila melanogaster"), c("Hs", "Mm", "Rn", "Dm"))
    if(what == "scientific"){
        return(names[names[,2] == organism, 3])
    }else{
        return(names[names[,1] == organism, 3])
    }
}

#getFullIDName <- function(ID){
#    switch(toupper(ID),
#           LL = return("LocusLink id"),
#           UG = return("UniGene id"),
#           GO = return("Gene Ontology id"),
#           PMID = return("PubMid id"),
#           stop("Unknown id name"))
#}

# Just do not what to throw the code away
.keepAround <- function(){
    # Read data with mappings between LL and GG
    go <- read.table(paste(url, getExten("go"), sep = ""), sep = "\t",
                        header = FALSE, as.is = TRUE)
    ll2GO <- mergeRowByKey(go, keyCol = 1)
    saveColSepData(ll2GO, pkgName, pkgPath, "LL2GO")
    copySubstitute(file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "#REPLACEME#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "LL2GO.Rd", sep = "")),
                   list(REPLACEME = "LL2GO", KEY = "LocusLink id",
                        VALUE = "Gene Ontology id"), "#")
    go2LL <- mergeRowByKey(go, keyCol = 2)
    saveColSepData(go2LL, pkgName, pkgPath, "GO2LL")
    copySubstitute(file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "#REPLACEME#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "GO2LL.Rd", sep = "")),
                   list(REPLACEME = "GO2LL", VALUE = "LocusLink id",
                        KEY = "Gene Ontology id"), "#")
    # Read data with mappings between LL and RefSeq
    ref <-  read.table(paste(url, getExten("ref"), sep = ""), sep = "\t",
                        header = FALSE, as.is = TRUE)
    ll2Ref <- mergeRowByKey(cbind(ref[,1], gsub("(^.*)\\..*", "\\1",
                                                   ref[,2])), keyCol = 1)
    saveColSepData(ll2Ref, pkgName, pkgPath, "LL2REF")
    copySubstitute(file.path(pkgPath, pkgName, "man",
                              paste(pkgName, "#REPLACEME#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "LL2REF.Rd", sep = "")),
                   list(REPLACEME = "LL2REF", KEY = "LocusLink id",
                        VALUE = "RefSeq id"), "#")
    ref2LL <- mergeRowByKey(cbind(ref[,1], gsub("(^.*)\\..*", "\\1",
                                                   ref[,2])), keyCol = 2)
    saveColSepData(ref2LL, pkgName, pkgPath, "REF2LL")
    copySubstitute(file.path(pkgPath, pkgName, "man",
                              paste(pkgName, "#REPLACEME#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "REF2LL.Rd", sep = "")),
                   list(REPLACEME = "REF2LL", VALUE = "LocusLink id",
                        KEY = "RefSeq id"), "#")
    # Read data with mappings between LL and UG
    ug <-  read.table(paste(url, getExten("ug"), sep = ""), sep = "\t",
                        header = FALSE, as.is = TRUE)
    ll2UG <- mergeRowByKey(ug, keyCol = 1)
    saveColSepData(ll2UG, pkgName, pkgPath, "LL2UG")
    copySubstitute(file.path(pkgPath, pkgName, "man",
                              paste(pkgName, "#REPLACEME#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "LL2UG.Rd", sep = "")),
                   list(REPLACEME = "LL2UG", KEY = "LocusLink id",
                        VALUE = "UniGene id"), "#")
    ug2LL <- mergeRowByKey(ug, keyCol = 2)
    saveColSepData(ug2LL, pkgName, pkgPath, "UG2LL")
    copySubstitute(file.path(pkgPath, pkgName, "man",
                              paste(pkgName, "#REPLACEME#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "UG2LL.Rd", sep = "")),
                   list(REPLACEME = "UG2LL", VALUE = "LocusLink id",
                        KEY = "UniGene id"), "#")
    # Write man pages and so on
    writeDescription(pkgName, pkgPath, version, author)
    writeFun (pkgPath, pkgName, organism = "")
    writeMan4Fun(pkgName, pkgPath)
    # Write .First.lib
    writeZZZ(pkgPath, pkgName)
    # Write the quality control data
    # Quality control
    getDPStats("", pkgName, pkgPath)
    writeMan4QC(pkgName, pkgPath)
    file.remove(file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "#REPLACEME#.Rd", sep = "")))
}

getExten <- function(what){

    switch(toupper(what),
           "GO" = return("loc2go"),
           "REF" = return("loc2ref"),
           "UG" = return("loc2UG"),
           "LL" = return("LL_tmpl.gz"),
           "ACC" = return("loc2acc"),
           stop("Source type not supported"))
}

saveData2Env <- function(data, fun = splitEntry, pkgName, pkgPath, envName){
#    data <- data[!is.na(data[,2]),]
    data <- data[data[,2]!="NA",]
    env <- new.env(hash = TRUE, parent = NULL)
    if(nrow(data) != 0){
        multiassign(data[,1], lapply(data[,2], fun), env)
    }
    assign(paste(pkgName, envName, sep = ""), env)
    save(list = paste(pkgName, envName, sep = ""),
         file = file.path(pkgPath, pkgName, "data",
         paste(pkgName, envName, ".rda", sep = "")))
}


getReverseMapping <- function(data, sep = ";"){
    if(nrow(data) == 0){
        return(c(data[2], data[1]))
    }
    seped <- sapply(data[,2], splitEntry, sep = sep)
    reps <- sapply(seped, length)
    data <- cbind(rep(data[,1], reps), unlist(seped, use.names = FALSE))
    data <- mergeRowByKey(data, keyCol = 2)
    return(data)
}

reverseMap4GO <- function(data, sep = ";", type = c("ll2GO", "GO2LL")){
    type <- match.arg(type)

    appendEC <- function(entry){
        ids <- unlist(strsplit(entry[2], ";"))
        if(type == "ll2GO"){
            return(paste(paste(ids, gsub("GO:[0-9]*", "", entry[1]),
                           sep = ""), sep = "", collapse = ";"))
        }else{
            return(paste(paste(ids, gsub("[0-9]*", "", entry[1]),
                           sep = ""), sep = "", collapse = ";"))
        }
    }
    if(nrow(data) == 0){
        if(length(data) != 2){
            return(data)
        }else{
            data <- cbind(data[2], data[1])
        }
    }
    # Get reverse mappings first
    reversed <- getReverseMapping(data, sep)
    # Append "@evidence code" to the end of each mapped values
    reversed[,2] <- apply(reversed, 1, appendEC)

    # Remove evidence code from GO ids
    reversed[,1] <- gsub("@.*", "", reversed[,1])

    return(mergeRowByKey(reversed, keyCol = 1))
}


getLL2ACC <- function(url =
                      paste("ftp://ftp.ncbi.nih.gov/refseq/LocusLink/",
                            getExten("acc"), sep = ""), organism = "human"){
    conn <- url(url)
    ll2Acc <- read.delim(url, sep = "\t", header = FALSE,
                                   strip.white = TRUE, as.is = TRUE)
    close(conn)

    orgNameCode <- getOrgNameNCode()

    ll2Acc <- ll2Acc[ll2Acc[,6] == names(orgNameCode[
                                     orgNameCode == getOrgName(organism,
                                     "scientific")]), 1:2]
    # The source data contain some "none"s, drop them
    ll2Acc <- ll2Acc[ll2Acc[,2] != "none", ]
    # Drop the version numbers
    ll2Acc[,2] <- gsub("\\..*$", "", ll2Acc[, 2])

    colnames(ll2Acc) <- c("LOCUSID", "ACCNUM")

    return(mergeRowByKey(ll2Acc))
}
