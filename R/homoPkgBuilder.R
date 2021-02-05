# This function builds a data package that maps internal HomoloGene
# ids of an organism to LocusLink ids, UniGene ids, percent identity
# of the alignment, and type of similarities of organisms of all
# pairwise best matches based on data from
# "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/hmlg.ftp"
# pkgName - a character string for the name of the data package to be
# created.
# pkgPath - a character string for the path to which the data package
# to be created will be stored.
# version - a character string for the version number of the data
# package to be created.
# author - a list of two elements with the firest one being author - a
# character string for the name of the author of the data package and
# second one being maintainer - a character string for the email of
# the author.
#

homoPkgBuilder <- function(suffix = "homology", pkgPath, version, author,
                 url = getSrcUrl("HG")){

    makeSrcInfo()
    homoMapping <- procHomoData(url)
    # The source data only gives one way relation. Include the reverse
    # relations too
    homoMapping <- as.matrix(rbind(homoMapping[,
                           c(1, 4, 5, 6, 2, 7, 8, 9, 10, 3)],
                           homoMapping[, c(2, 7, 8, 9, 1, 4, 5, 6, 10, 3)]))
    # Some relations contain genes with unknown LLids. Keep them under 0
    #homoMapping[is.na(homoMapping[,2]), 2] <- "0"
    tempList <- getRepList("HG", list(hg = HG()))
    for(i in unique(homoMapping[, 1])){
        print(paste("Processing data for", mapOrgs(i, "code")))
        pkgName <- paste(tolower(getShortSciName(mapOrgs(i, "code"))),
                              suffix, sep = "")
        createEmptyDPkg(pkgName, pkgPath)
        tempHomo <- homoMapping[homoMapping[,1] == i, ]
        tempList[["PKGNAME"]] <- pkgName
        tempList[["ORGCODE"]] <- i
        tempList[["ORGANISM"]] <- mapOrgs(i, "code")
        mapPS(tempHomo, pkgName, pkgPath, tempList)
        writeDescription(pkgName, pkgPath, version, author)
        writeFun (pkgPath, pkgName, organism = "")
        writeZZZ(pkgPath, pkgName)
        # Write the quality control data
        getDPStats("", pkgName, pkgPath)
        copySubstitute(file.path(.path.package("AnnBuilder"),
                        "templates", "PKGNAMEQC.Rd"),
                        file.path(pkgPath, pkgName, "man",
                        paste(pkgName, "QC.Rd", sep = "")),
                        list(PKGNAME = pkgName), "#")
        copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "homology.Rd"),
                       file.path(pkgPath, pkgName,
                                 "man", paste(pkgName, ".Rd", sep = "")),
                                 tempList, "#")
        makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
    }
}

procHomoData <- function(url = getSrcUrl("HG")){
    # scan() does not work for some reason
    homoData <- as.matrix(read.table(url, header = FALSE, sep = "|",
			stringsAsFactors=FALSE,
                           quote = "", as.is = TRUE, strip.white = TRUE,
                           comment.char = ""))
    # remove leading LL from LL ids
    homoData <- matrix(gsub("^LL\\.", "", homoData),
                       ncol = ncol(homoData), byrow = FALSE)
    # make "" NA
    homoData[homoData == ""] <- NA
    # remove leading white spaces
    homoData <- matrix(gsub(" *", "", homoData),
                       ncol = ncol(homoData), byrow = FALSE)
    return(homoData)
}

mapPS <- function(homoMappings, pkgName, pkgPath, tempList){
    #homoMapping <- procHomoData(url)
    # The source data only gives one way relation. Include the reverse
    # relations too
    #homoMapping <- as.matrix(rbind(homoMapping[,
    #                       c(1, 4, 5, 6, 2, 7, 8, 9, 10, 3)],
    #                       homoMapping[, c(2, 7, 8, 9, 1, 4, 5, 6, 10, 3)]))
    # Some relations contain genes with unknown LLids. Keep them under 0
    #homoMapping[is.na(homoMapping[,2]), 2] <- "0"
    #tempList <- getRepList("HG", list(hg = HG()))
    # Get the vector of HGIDs for each organism
    #print("Processing HGIDs for each organism")
    #for(i in unique(homoMapping[, 1])){
        #tempHomo <- homoMappings[homoMappings[,1] == i, ]
        tempVect <- unique(as.vector(homoMappings[, 3]))
        tempVect <- tempVect[!is.na(tempVect)]
        assign(paste(pkgName, "HGID", sep = ""), tempVect)
        save(list = paste(pkgName, "HGID", sep = ""),
             file = file.path(pkgPath, pkgName, "data",
             paste(pkgName, "HGID", ".rda", sep = "")))
        copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "homologyHGID.Rd"),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "HGID.Rd",
                                        sep = "")), tempList, "#")
        #tempList[["ORGCODE"]] <- i
        #tempList[["ORGANISM"]] <- mapOrgs(i, "code")

        #print("Creating lists of homoData objects")
        homodata <- HomoData2List(homoMappings[, c(3, 5:ncol(homoMappings))],
                                  what = "old")
        saveList(homodata, pkgName, pkgPath, "DATA")
        copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "homologyDATA.Rd"),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "DATA.Rd",
                                        sep = "")), tempList, "#")
        #print("Mapping LL ids to HGIDs")
        ll2HGID <- cbind(homoMappings[, 2], paste(homoMappings[, 3],
                                          homoMappings[, 1], sep = "@"))
        ll2HGID <- ll2HGID[!is.na(ll2HGID[,1]), ]
        ll2HGID <- mergeRowByKey(ll2HGID)
        saveMat(ll2HGID, pkgName, pkgPath, "LL2HGID",
                fun = function(x) twoStepSplit(x, asNumeric = TRUE))
        copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "homologyLL2HGID.Rd"),
                           file.path(pkgPath, pkgName, "man",
                                     paste(pkgName, "LL2HGID.Rd", sep = "")),
                                     tempList, "#")
        #print("Mapping HGIDs to LL")
        HGID2LL <- cbind(homoMappings[, 3], paste(homoMappings[, 2],
                                          homoMappings[, 1], sep = "@"))
        HGID2LL <- mergeRowByKey(HGID2LL)
        saveMat(HGID2LL, pkgName, pkgPath, "HGID2LL",
            fun = function(x) twoStepSplit(x, asNumeric = TRUE))
        #print("Mapping Acc number to HGIDs")
        homoMappings <- homoMappings[!is.na(homoMappings[,4]),]
        acc2HGID <- cbind(homoMappings[,4], paste(homoMappings[, 3],
                                          homoMappings[, 1], sep = "@"))
        acc2HGID[,1] <- gsub("\\..*$", "", acc2HGID[,1])
        acc2HGID <- acc2HGID[!is.na(acc2HGID[,1]), ]
        acc2HGID <- mergeRowByKey(acc2HGID)
        saveMat(acc2HGID, pkgName, pkgPath, "ACC2HGID",
            fun = function(x) twoStepSplit(x, asNumeric = TRUE))
        #print("Mapping HGIDs to Acc numbers")
        HGID2ACC <- cbind(homoMappings[,3],
                      paste(gsub("\\..*", "", homoMappings[, 4]),
                            homoMappings[, 1], sep = "@"))
        HGID2ACC <- mergeRowByKey(HGID2ACC)
        saveMat(HGID2ACC, pkgName, pkgPath, "HGID2ACC", fun = twoStepSplit)

        saveOrgNameNCode(pkgName, pkgPath, tempList)
        #copySubstitute(file.path(.path.package("AnnBuilder"),
        #                         "templates", "homologyORGCODE.Rd"),
        #               file.path(pkgPath, pkgName,
        #                         "man", "homology.Rd"), tempList, "#")
    #}
}

HomoData2List <- function(data, what = "old"){
    data <- data[!is.na(data[,1]),]
    temp <- split.data.frame(data, factor(data[,1]))
    return(sapply(temp, getHomoDList, what = what))
}

getHomoDList <- function(data, what = "old"){
    temp <- apply(data, 1, getHomoData, what = what)
    names(temp) <- data[,2]
    return(temp)
}

getHomoData <- function(entries, what = "old", objOK = FALSE){
    #entries[entries == "NA"] <- NA
    if(what == "old"){
        if(length(entries) < 7){
            stop("Incorrect argument length. Must be a vector of 7")
        }
        if(!is.na(entries[7]) && any(entries[7] == c("B", "b", "m"))){
            if(objOK){
                return(new("homoData",
                           homoOrg =
                           tolower(getShortSciName(mapOrgs(entries[2]))),
                           #omoLL = as.numeric(entries[3]),
                           homoType = as.character(entries[7]),
                           homoPS = as.numeric(entries[6]),
                           #omoACC = as.character(entries[5]),
                           homoHGID = as.numeric(entries[4]),
                           homoURL = as.character(NA)))
            }else{
                return(list(homoOrg =
                            tolower(getShortSciName(mapOrgs(entries[2]))),
                            #omoLL = as.numeric(entries[3]),
                            homoType = as.character(entries[7]),
                            homoPS = as.numeric(entries[6]),
                            #omoACC = as.character(entries[5]),
                            homoHGID = as.numeric(entries[4]),
                            homoURL = as.character(NA)))
            }
        }else{
            if(objOK){
                return(new("homoData",
                           homoOrg =
                           as.character(tolower(getShortSciName(mapOrgs(entries[2])))),
                           #omoLL = as.numeric(entries[3]),
                           homoType = as.character(entries[7]),
                           #omoACC = as.character(entries[5]),
                           homoHGID = as.numeric(entries[4]),
                           homoPS = as.numeric(NA),
                           homoURL = as.character(entries[6])))
            }else{
                return(list(homoOrg =
                            as.character(tolower(getShortSciName(mapOrgs(entries[2])))),
                           #omoLL = as.numeric(entries[3]),
                           homoType = as.character(entries[7]),
                           #omoACC = as.character(entries[5]),
                           homoHGID = as.numeric(entries[4]),
                           homoPS = as.numeric(NA),
                           homoURL = as.character(entries[6])))
            }
        }
    }else{
        if(length(entries) < 11){
            stop("Incorrect argument length. Must be a vector of 10")
        }
        entries <- gsub(" ", "", entries)
        if(objOK){
            return(.newHomoData(as.character(tolower(getShortSciName(mapOrgs(entries[3])))),
                            #as.numeric(entries[2]),
                            homoBest = as.logical(entries[4]),
                            homoNucChange = as.numeric(entries[5]),
                            homoNucChangeJc = as.numeric(entries[6]),
                            homoProtChange = as.numeric(entries[7]),
                            homoKA = as.numeric(entries[8]),
                            homoKS = as.numeric(entries[9]),
                            homoKNR = as.numeric(entries[10]),
                            homoKNC = as.numeric(entries[11])))
        }else{
            return(list(homoOrg =
                        as.character(tolower(getShortSciName(mapOrgs(entries[3])))),
                        #homoLL = as.numeric(entries[2]),
                        homoBest = as.logical(entries[4]),
                        homoNucChange = as.numeric(entries[5]),
                        homoNucChangeJc = as.numeric(entries[6]),
                        homoProtChange = as.numeric(entries[7]),
                        homoKA = as.numeric(entries[8]),
                        homoKS = as.numeric(entries[9]),
                        homoKNR = as.numeric(entries[10]),
                        homoKNC = as.numeric(entries[11])))
        }
    }
}

saveOrgNameNCode <- function(pkgName, pkgPath, tepList){
    tepList[["PKGNAME"]] <- pkgName
    temp <- getOrgNameNCode()
    #temp <- cbind(temp, names(temp), tolower(sapply(temp, getShortSciName)))
    temp <- data.frame(cbind(as.matrix(temp), as.matrix(names(temp)), as.matrix(tolower(sapply(temp, getShortSciName)))))
    colnames(temp) <- c("species_name","tax_id","tla")
    rownames(temp) <- c(1:nrow(temp))
    #codeNName <- names(temp)
    #names(codeNName) <- temp
    assign(paste(pkgName, "ORGCODE", sep = ""), temp)
    save(list = paste(pkgName, "ORGCODE", sep = ""),
    	file =  file.path(pkgPath, pkgName, "data",
    	paste(pkgName, "ORGCODE.rda", sep = "")))
    copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "homologyORGCODE.Rd"),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, "ORGCODE.Rd", sep = "")),
                                     tepList, "#")
}

getLL2IntID <- function(homoData, organism = ""){
    switch(tolower(organism),
           human = org <- "9606",
           mouse = org <- "10090",
           rat = org <- "10116",
           org <- NULL)
    ll2ID <- homoData[, c(4, 5, 1)]
    ll2ID <- rbind(ll2ID, homoData[,c(7, 8, 2)])
    if(!is.null(org)){
        ll2ID <- ll2ID[ll2ID[,3] == org, ]
    }
    ll2ID <- ll2ID[, 1:2]
    ll2ID <- ll2ID[ll2ID[, 1] != "NA", ]

    if(nrow(ll2ID) > length(unique(ll2ID[,1]))){
        ll2ID <- mergeRowByKey(ll2ID)
    }
    colnames(ll2ID) <- c("ENTREZID", "HGID")

    return(ll2ID)
}

writeHomoXMLData <- function(pkgName = "homology", pkgPath, version, author,
    url = "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build39.2/homologene.xml.gz"){

    homoXML <- loadFromUrl(url, tempdir())
    homodata <- homoXMLParser(homoXML)
    unlink(homoXML)

    makeSrcInfo()
    tempList <- getRepList("HG", list(hg = HG()))
    createEmptyDPkg(pkgName, pkgPath)
    writeHGID2Caption(pkgName, pkgPath, homodata[["hgid2cap"]])
    #writeOrgLL(pkgName, pkgPath, homodata[["ll2Org"]])
    writeHGID2LL(pkgName, pkgPath, homodata[["hgid2gi"]])
    writeHomoData(pkgName, pkgPath, homodata[["homodata"]])

    writeDescription(pkgName, pkgPath, version, author)
    writeFun (pkgPath, pkgName, organism = "")
    writeZZZ(pkgPath, pkgName)
    tempList[["PKGNAMW"]] <- pkgName
    saveOrgNameNCode(pkgName, pkgPath, tempList)
    # Write the quality control data
    getDPStats("", pkgName, pkgPath)
    copySubstitute(file.path(.path.package("AnnBuilder"),
                        "templates", "PKGNAMEQC.Rd"),
                        file.path(pkgPath, pkgName, "man",
                        paste(pkgName, "QC.Rd", sep = "")),
                        list(PKGNAME = pkgName), "#")
    writeMan4QC(pkgName, pkgPath)

    return(invisible(NA))
}

writeHGID2Caption <- function(pkgName, pkgPath, hgid2Cap){
    saveMat(hgid2Cap, pkgName, pkgPath, "HGID2CAPTION")
    saveMat(getReverseMapping(hgid2Cap), pkgName, pkgPath, "CAPTION2HGID")
    copySubstitute(file.path(.path.package("AnnBuilder"),
                             "templates", "homologyHGID2CAPTION.Rd"),
                   file.path(pkgPath, pkgName,
                             "man", "homologyHGID2CAPTION.Rd"),
                   getRepList("HG", list(hg = HG())), "#")
}

#writeOrgLL <- function(pkgName, pkgPath, ll2Org){
#    tempList <- getRepList("HG", list(hg = HG()))
#    print("Processing LLIDs for each organism")
#    for(i in unique(ll2Org[, 2])){
#        tempVect <- unique(as.vector(ll2Org[ll2Org[, 2] == i, 1]))
#        tempVect <- tempVect[!is.na(tempVect)]
#        assign(paste(pkgName, i, sep = ""), tempVect)
#        save(list = paste(pkgName, i, sep = ""),
#             file = file.path(pkgPath, pkgName, "data",
#             paste(pkgName, i, ".rda", sep = "")))
#        tempList[["ORGCODE"]] <- i
#        tempList[["ORGANISM"]] <- mapOrgs(i, "code")
#        copySubstitute(file.path(.path.package("AnnBuilder"),
#                                 "templates", "homology#ORGCODE#.Rd"),
#                           file.path(pkgPath, pkgName,
#                           "man", paste(pkgName, i, ".Rd",
#                                        sep = "")), tempList, "#")
#    }
#}

writeHGID2LL <- function(pkgName, pkgPath, hgid2LL){
    hgid2LL <- hgid2LL[!is.na(hgid2LL[,1]), ]
    saveMat(hgid2LL, pkgName, pkgPath, "HGID2LL", fun = twoStepSplit)
    ll2HGID <- sapply(hgid2LL[,2], splitEntry)
    reps <- sapply(ll2HGID, length)
    ll2HGID <- cbind(rep(hgid2LL[,1], reps),
                          unlist(ll2HGID, use.names = FALSE))
    ll2HGID[,2] <- paste(ll2HGID[1,], "@",
                         gsub(".*(@.*$)", "\\1", ll2HGID[1,]))
    ll2HGID <- mergeRowByKey(ll2HGID)
    saveMat(ll2HGID, pkgName, pkgPath, "LL2HGID", fun = twoStepSplit)
    copySubstitute(file.path(.path.package("AnnBuilder"),
                             "templates", "homologyLL2HGID.Rd"),
                   file.path(pkgPath, pkgName,
                             "man", "homologyLL2HGID.Rd"),
                   getRepList("HG", list(hg = HG())), "#")
}

writeHomoData <- function(pkgName, pkgPath, homoFile){

    #getHomoData <- function(id){
    #    fun <- function(entry){
    #        entry <- gsub(" *", "", as.vector(entry))
    #        return(homoData(as.character(entry[2]),
    #                  as.numeric(entry[1]), as.logical(entry[3]),
    #                  as.numeric(entry[4]), as.numeric(entry[5]),
    #                  as.numeric(entry[6]), as.numeric(entry[7]),
    #                  as.numeric(entry[8]), as.numeric(entry[9]),
    #                  as.numeric(entry[10])))
    #    }
    #    data <- homodata[homodata[,1] == id, c(2, 4, 5:12)]
    #    tempList <- apply(data, 1, fun)
    #    names(tempList) <- gsub(" *", "", as.vector(data[,1]))
    #    return(tempList)
    # }

    homodata <- read.table(homoFile, sep = "\t", header = FALSE,
			stringsAsFactors=FALSE,
                           as.is = TRUE, strip.white = TRUE)
    # homodata has one way homology mappings. Add the mappings that
    # are reciprocal.
    recip <- homodata[homodata[,5] == "true", c(2, 1, 4, 3, 5:12)]
    colnames(recip) <- colnames(homodata)
    homodata <- rbind(homodata, recip)
    tempList <- getRepList("HG", list(hg = HG()))
    for(i in unique(homodata[,3])){
        print(paste("Processing data for", mapOrgs(as.character(i))))
        tempHomoData <- HomoData2List(homodata[homodata[,3] == i,
                                   c(1:2, 4:ncol(homodata))], what = "xml")
        saveList(tempHomoData, pkgName, pkgPath, gsub(" ", "", i))
        tempList[["ORGCODE"]] <- gsub(" ", "", i)
        tempList[["ORGANISM"]] <- mapOrgs(gsub(" ", "", i), "code")
        copySubstitute(file.path(.path.package("AnnBuilder"),
                                 "templates", "homology#ORGCODE#.Rd"),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, gsub(" ", "", i), ".Rd",
                                        sep = "")), tempList, "#")
    }

    #homoList <- HomoData2List(homodata[, c(1:2, 4:ncol(homodata))], "xml")

    #saveList(homoList, pkgName, pkgPath, "homoDATA")
    #copySubstitute(file.path(.path.package("AnnBuilder"),
    #                         "templates", "homologyHOMODATA.Rd"),
    #               file.path(pkgPath, pkgName,
    #                         "man", "homologyHOMODATA.Rd"),
    #               getRepList("HG", list(hg = HG())), "#")
    unlink(homoFile)
    return(invisible())
}

homoXMLParser <- function(fileName) {
    text <- NA
    hgid <- NA
    caption <- NA
    gi <- NA
    org <- NA
    protGi <- NA
    id1 <- NA
    id2 <- NA
    nucChange <- NA
    nucChangeJc <- NA
    protChange <- NA
    ka <- NA
    ks <- NA
    knr <- NA
    knc <- NA
    best <- NA


    hgid2Caption <- matrix(nrow = 0, ncol = 2)
    hgid2Gi <- matrix(nrow = 0, ncol = 2)
    gene <- matrix(nrow = 0, ncol = 3)
    homoFile <- tempfile()
    ll2Org <- matrix(nrow = 0, ncol = 2)

    reSetValues <- function(){
        text <<- NA
        hgid <<- NA
        caption <<- NA
        gi <<- NA
        org <<- NA
        protGi <<- NA
        gene <<- matrix(nrow = 0, ncol = 3)
        reSetGene()
    }

    reSetGene <- function(){
        id1 <<- NA
        id2 <<- NA
        nucChange <<- NA
        nucChangeJc <<- NA
        protChange <<- NA
        ka <<- NA
        ks <<- NA
        knr <<- NA
        knc <<- NA
        best <<- NA
    }

    writeEntry <- function(){
        if(!is.na(hgid)){
            hgid2Caption <<- rbind(hgid2Caption, c(hgid, caption))
            hgid2Gi <<- rbind(hgid2Gi, c(hgid, paste(gene[, 1], gene[,2],
                                             sep = "@", collapse = ";")))
            ll2Org <<- rbind(ll2Org, gene[, 1:2])
        }
        reSetValues()
    }

    writeGene <- function(){
        gene <<- rbind(gene, c(gi, org, protGi))
        gi <<- NA
        org <<- NA
        protGi <<- NA
    }

    writeDist <- function(){
        gi1 <- gene[gene[,3] == id1, 1]
        gi2 <- gene[gene[,3] == id2, 1]
        homodata <- paste(gi1, gi2, gene[gene[,1] == gi1, 2],
                    gene[gene[,1] == gi2, 2], best, nucChange,
                    nucChangeJc, protChange, ka, ks, knr, knc,
                    sep = "\t")
        write(homodata, file = homoFile, append = TRUE)

        #if(!gi %in% names(homodata)){
        #    tempList <- list()
        #    tempList[[gene[gene[,3] == id2, 1]]] <-
        #                            getHomoData(gene[gene[,3] == id2, 1])
        #    homodata[[gi]] <<- tempList
        #}else{
        #    if(!gene[gene[,3] == id2, 1] %in%
        #       unlist(sapply(homodata[[gi]], FUN = function(x) homoLL(x)))){
        #        homodata[[gi]][[gene[gene[,3] == id2, 1]]] <<-
        #                             getHomoData(gene[gene[,3] == id2, 1])
        #    }
        #}
        reSetGene()
    }

    # Event handlers for the parser
    getHandler <- function(){
        # Gets the value for desired element
        startElement <- function(name, attrs){
            switch(name,
                   "HomoloGeneEntry" = writeEntry(),
                   "Stats_recip-best" = best <<- attrs
                  )
        }
        # Write the data when entries for a GO term ends
        endElement <- function(name,...){
            switch(name,
                   "HomoloGeneEntry_hg-id" = hgid <<- text,
                   "HomoloGeneEntry_caption" = caption <<- text,
                   "Gene_geneid" = gi <<- text,
                   "Gene_taxid" = org <<- text,
                   "Gene_prot-gi" = protGi <<- text,
                   "Gene" = writeGene(),
                   "Stats_gi1" = id1 <<- text,
                   "Stats_gi2" = id2 <<- text,
                   "Stats_nuc-change" = nucChange <<- text,
                   "Stats_nuc-change-jc" = nucChangeJc <<- text,
                   "Stats_prot-change" = protChange <<- text,
                   "Stats_ka" = ka <<- text,
                   "Stats_ks" = ks <<- text,
                   "Stats_knr" = knr <<- text,
                   "Stats_knc" = knc <<- text,
                   "Stats" = writeDist()
                   )
        }
        text <- function(x){
            text <<- x
        }

        list(startElement = startElement, endElement = endElement,
             text = text)
    }
    xmlEventParse(fileName, handlers = getHandler())

    return(list(hgid2cap = hgid2Caption, hgid2gi = hgid2Gi,
                ll2Org = ll2Org[!duplicated(ll2Org),],
                homodata = homoFile))
}


# Keep the new definition for now
.newHomoData <- function(homoOrg, homoLL, homoBest, homoNucChange,
                         homoNucChangeJc, homoProtChange, homoKA,
                         homoKS, homoKNR, homoKNC){
  setClass("homoData", representation(homoOrg = "character",
                                    homoLL = "numeric",
                                    homoBest = "logical",
                                    homoNucChange = "numeric",
                                    homoNucChangeJc = "numeric",
                                    homoProtChange = "numeric",
                                    homoKA = "numeric",
                                    homoKS = "numeric",
                                    homoKNR = "numeric",
                                    homoKNC = "numeric"))

# Set the get methods

    setGeneric("homoOrg",
               function(object) standardGeneric("homoOrg"))

setMethod("homoOrg", "homoData",
          function(object) object@homoOrg)


    setGeneric("homoLL",
               function(object) standardGeneric("homoLL"))

setMethod("homoLL", "homoData",
          function(object) object@homoLL)


    setGeneric("homoBest",
               function(object) standardGeneric("homoBest"))

setMethod("homoBest", "homoData",
          function(object) object@homoBest)


    setGeneric("homoNucChange",
               function(object) standardGeneric("homoNucChange"))

setMethod("homoNucChange", "homoData",
          function(object) object@homoNucChange)


    setGeneric("homoNucChangeJc",
               function(object) standardGeneric("homoNucChangeJc"))

setMethod("homoNucChangeJc", "homoData",
          function(object) object@homoNucChangeJc)


    setGeneric("homoProtChange",
               function(object) standardGeneric("homoProtChange"))

setMethod("homoProtChange", "homoData",
          function(object) object@homoProtChange)


    setGeneric("homoKA",
               function(object) standardGeneric("homoKA"))

setMethod("homoKA", "homoData",
          function(object) object@homoKA)


    setGeneric("homoKS",
               function(object) standardGeneric("homoKS"))

setMethod("homoKS", "homoData",
          function(object) object@homoKS)


    setGeneric("homoKNR",
               function(object) standardGeneric("homoKNR"))

setMethod("homoKNR", "homoData",
          function(object) object@homoKNR)


    setGeneric("homoKNC",
               function(object) standardGeneric("homoKNC"))

setMethod("homoKNC", "homoData",
          function(object) object@homoKNC)

setMethod("show", "homoData",
          function(object) {
              if(length(homoOrg(object)) > 0 && !is.na(homoOrg(object))){
                  cat(paste("homoOrg:", homoOrg(object)))
              }
              if(length(homoLL(object)) > 0 && !is.na(homoLL(object))){
                  cat(paste("\nhomoLL:", homoLL(object)))
              }
              if(length(homoBest(object)) > 0 && !is.na(homoBest(object))){
                  cat(paste("\nhomoBest:", homoBest(object)))
              }
              if(length(homoNucChange(object)) > 0 &&
                                    !is.na(homoNucChange(object))){
                  cat(paste("\nhomoNucChange:", homoNucChange(object)))
              }
              if(length(homoNucChangeJc(object)) > 0 &&
                                  !is.na(homoNucChangeJc(object))){
                  cat(paste("\nhomoNucChangeJc:", homoNucChangeJc(object)))
              }
              if(length(homoProtChange(object)) > 0 &&
                                  !is.na(homoProtChange(object))){
                  cat(paste("\nhomoProtChange:", homoProtChange(object)))
              }
              if(length(homoKA(object)) > 0 && !is.na(homoKA(object))){
                  cat(paste("\nhomoKA:", homoKA(object)))
              }
              if(length(homoKS(object)) > 0 && !is.na(homoKS(object))){
                  cat(paste("\nhomoKS:", homoKS(object)))
              }
              if(length(homoKNR(object)) > 0 && !is.na(homoKNR(object))){
                  cat(paste("\nhomoKNR:", homoKNR(object)))
              }
              if(length(homoKNC(object)) > 0 && !is.na(homoKNC(object))){
                  cat(paste("\nhomoKNC:", homoKNC(object)))
              }
              cat("\n")
             
              
            })


  return(new("homoData", homoOrg = homoOrg,
                   homoLL = homoLL, homoBest = homoBest,
                   homoNucChange = homoNucChange,
                         homoNucChangeJc = homoNucChangeJc,
                   homoProtChange = homoProtChange, homoKA = homoKA,
                   homoKS = homoKS, homoKNR = homoKNR, homoKNC = homoKNC))
}



