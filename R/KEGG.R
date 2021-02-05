## Constructor for KEGG objects
KEGG <- function(srcUrl = "ftp://ftp.genome.ad.jp/pub/kegg/pathway/organisms",
                 organism = "human", parser = "", baseFile = "",
                 built, fromWeb = TRUE,
                 kegggenomeUrl = "ftp://ftp.genome.ad.jp/pub/kegg/genes/genome"){

    KEGGNames <- parseKEGGGenome(url=kegggenomeUrl)
    assignInNamespace("KEGGNames", KEGGNames, "AnnBuilder")
##           pos = match("package:AnnBuilder", search()))

    new("KEGG", srcUrl = srcUrl, organism = organism, parser = parser,
        baseFile = baseFile, built = ifelse(missing(built),
                             getSrcBuilt("KEGG", organism), built),
        fromWeb = fromWeb)

}


## An placeholder so we can assign into the namespace for this binding
KEGGNames <- NULL

getKEGGIDNName <- function(object, exten = "/../map_title.tab" ){
    fullUri <- paste(srcUrl(object), exten, sep = "")
    if(fromWeb(object)){
        con <- url(fullUri, "r")
    }else{
        con <- file(fullUri)
    }
    temp <- matrix(unlist(strsplit(readLines(con), "\t")),
                                      ncol = 2, byrow = TRUE)
    close(con)
    idNName <- temp[,2]
    names(idNName) <- temp[,1]
    return(idNName)
}

getKEGGOrgName <- function(name){
    options(show.error.messages = FALSE)
##     orgName <- get("KEGGNames", pos = match("package:AnnBuilder",
##                     search()))[tolower(KEGGNames[,2]) == tolower(name), 1]
    orgName <- AnnBuilder:::KEGGNames[tolower(KEGGNames[,2]) == tolower(name), 1]
    options(show.error.messages = TRUE)
    if(inherits(orgName, "try-error")){
        KEGGNames <- parseKEGGGenome()
        return(KEGGNames[tolower(KEGGNames[,2]) == tolower(name), 1])
    }else{
        return(orgName)
    }
}

parseKEGGGenome <- function(url){

    if(is.null(url)) {
        url="ftp://ftp.genome.ad.jp/pub/kegg/genes/genome"
    }
    
    keggURL <- url(url)
    options(show.error.messages = FALSE)
    temp <- try(readLines(keggURL))
    close(keggURL)
    options(show.error.messages = TRUE)
    if(inherits(temp, "try-error")){
        stop("Faild to obtain KEGG organism code")
    }
    temp <- temp[grep("ENTRY|DEFINITION|TAXONOMY", temp)]
    ## some KEGG entries are missing DEFINITION
    ## the following while() can fix this problem
    while( (3*(max(grep("DEFINITION", temp[(1:length(temp))%%3==2]))-1)+2) != (length(temp)-1) ) {
        temp <- c(temp[1:(3*(max(grep("DEFINITION", temp[(1:length(temp))%%3==2])))+1)],
                  "DEFINITION\t",
                  temp[(3*(max(grep("DEFINITION", temp[(1:length(temp))%%3==2])))+2):length(temp)]
                  )
    }
    
    temp <- matrix(temp, ncol = 3, byrow = TRUE)
    
    temp[,1] <- unlist(lapply(temp[,1], function(x)
                              unlist(strsplit(x, " +"))[2]))
    temp[,2] <- unlist(lapply(temp[,2], function(x)
                              paste(unlist(strsplit(x, " +"))[2:3],
                                    collapse = " ")))
    temp[,3] <- unlist(lapply(temp[,3], function(x)
                              unlist(strsplit(x, ":"))[2]))
    return(temp)
}

getLLPathMap <- function(srcUrl, idNName, organism, fromWeb = TRUE){
    llEC <- NULL
    llPathName <- NULL
    # Get ll to EC and pathway name mappings
    for(i in names(idNName)){

        temp <- mapll2EC(i, srcUrl = srcUrl, organism, fromWeb)
        if(!is.null(temp)){
            temp <- t(sapply(temp, parseEC))
            # Map ll to EC
            if(is.null(llEC)){
                llEC <- temp[!is.na(temp[,2]),]
            }else{
                llEC <- rbind(llEC, temp[!is.na(temp[,2]),] )
            }
            # Map ll to pathway name
            if(length(unique(temp[,1])) <= 1){
                tempPN <- matrix(c(unique(temp[,1]), i), ncol = 2)
            }else{
                tempPN <- cbind(unique(temp[, 1]), i)
            }
            if(is.null(llPathName)){
                llPathName <- tempPN
            }else{
                llPathName <- rbind(llPathName, tempPN)
            }
        }
    }
    if(is.null(ncol(llEC))){
        llEC <- matrix(NA, ncol = 2)
    }
    if(is.null(ncol(llPathName))){
        llPathName <- matrix(NA, ncol = 2)
    }
    colnames(llEC) <- c("ENTREZID", "ENZYME")
    colnames(llPathName) <- c("ENTREZID", "PATH")

    return(list(llec = mergeRowByKey(llEC),
                llpathname = mergeRowByKey(llPathName)))
}

mapll2PathID <- function(srcUrl, organism, exten = "gene_map.tab"){
    orgName <- getKEGGOrgName(organism)
    tempUrl <- paste(srcUrl, "/", orgName, "/", orgName, "_",
                     exten, sep = "")
    ll2Path <-as.matrix( read.table(tempUrl, header = FALSE, sep = "\t"))
    ll2Path <- cbind(ll2Path[, 1], gsub(" *", ";", ll2Path[,2]))
    colnames(ll2Path) <- c("ENTREZID", "PATH")

    return(ll2Path)
}

# Returns a matrix with Locus link ids as the first column and EC as
# second one with NA for ll not mapped to EC
mapll2EC <- function(id, srcUrl, organism, fromWeb, sep = "\t"){
        tempURL <- paste(srcUrl,"/", getKEGGOrgName(organism),
                     "/", getKEGGOrgName(organism), as.character(id),
                     ".gene", sep = "")
    if(fromWeb){
        con <- url(tempURL)
    }else{
        con <- file(tempURL)
    }

    on.exit(options(show.error.messages = TRUE))

    temp <- NULL
    #con <- url(tempURL)
    options(show.error.messages = FALSE)
    llNEC <- try(readLines(con))
    close(con)
    options(show.error.messages = TRUE)

    if(!inherits(llNEC, "try-error")){
        temp <- strsplit(llNEC, split = sep)
#        temp <- t(sapply(temp, parseEC))
#        temp <- temp[!is.na(temp[,2]),]
        return(temp)
    }else{
        # Turned off as a given organism may not have all the pathways 
        # cat("Failed to get data from URL:", tempURL, "\n")
        return(NULL)
    }
}

parseEC <- function(llNEC){
    if(length(grep(".*\\[EC:[0-9|\\.|-| ]*", llNEC[2])) != 0){
        EC <- strsplit(strsplit(llNEC[2], "\\[EC:")[[1]][2], "\\]")[[1]][1]
        # Collapse multiple EC number
        EC <- paste(unlist(strsplit(EC, " "), use.names = FALSE),
                    collapse = ";")

        return(c(llNEC[1], EC))
    }else{
        return(c(llNEC[1], NA))
    }
}

getKEGGGeneMap <- function(organism = "Homo sapiens"){
    idMap <- NULL
    for(i in organism){
        temp <- read.table(getKEGGFile("geneMap", getKEGGOrgName(i)),
			stringsAsFactors=FALSE,
                       sep = "\t", header = FALSE, as.is = TRUE)
        # append kegg orgnism id to the begining of each
        # KEGG pathway id and collapse multiple mappings
        temp <- cbind(temp[,1], unlist(lapply(paste(getKEGGOrgName(i),
                    temp[,2], sep = ""), gsub, pattern = " ",
                    replacement = paste(";", getKEGGOrgName(i), sep = ""))))
        idMap <- rbind(idMap, temp)
    }

    idMap <- idMap[idMap[,1] != "" & !is.na(idMap[,1]),]

    return(mergeRowByKey(idMap))
}













