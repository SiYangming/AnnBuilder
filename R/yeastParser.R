
yeastParser <- function(yeastData,
                         colKeep = c(1, 4, 5, 12)){
    domain <- read.delim(yeastData, header=FALSE, as.is=TRUE)
    domain <- domain[,colKeep]
    colnames(domain) <- c("PROBE", "METHOD", "ID", "INTERPROTMP")

    ##remove version #s from ID

    pasteParts = function(id,retId=""){
     for (i in 1:length(id)){
      retId = paste(retId,id[i],sep=".")
     }
      retId = sub(pattern="^\\.", replacement="", x=retId, perl=TRUE);
    }

    splitStrip = function(l){
     split = strsplit(l, "\\.")
     if(length(unlist(split[1]))>1){
      id = unlist(split[1])[1:length(unlist(split[1]))-1]
     }
     else{
      id = unlist(split[1])
     }
     id = toupper(pasteParts(id))
    }

#    domain$"ID" <- sapply(domain$"ID", function(x) unlist(strsplit(x, ".", fixed=TRUE))[1])
    domain$"ID" <- sapply(domain$"ID", splitStrip)

    domain$"INTERPROTMP" <- sapply(domain$"INTERPROTMP", function(x) unlist(strsplit(x, ".", fixed=TRUE))[1])

    ##indeies of methods
    isDomain <- list()
    isDomain[["PFAM"]]     <- domain$"METHOD" == "HMMPfam"
    isDomain[["SMART"]]    <- domain$"METHOD" == "HMMSmart"
    isDomain[["INTERPRO"]] <- domain$"INTERPROTMP" != "NULL"
    
    merged <- as.matrix(unique(domain$"PROBE"))
    colnames(merged) <- "PROBE"
    for(id in names(isDomain)){

        if(id != "INTERPRO"){
            idindex <- 3
        }else{
            idindex <- 4
        }

        tmpDomain <- domain[isDomain[[id]], c(1, idindex)]
        
        tmpList <- split(tmpDomain[, 2], tmpDomain[, 1])
        tmpList <- lapply(tmpList, unique)
        tmpList <- lapply(tmpList, function(x) paste(x, sep="", collapse=";"))
        tmp <- unlist(tmpList)
        tmp <- cbind(names(tmp), tmp)
        colnames(tmp) <- c("PROBE", id)
        merged <- merge(merged, tmp, by="PROBE", all.x=TRUE)

    }

    merged[merged==""] <- NA
    
    return(as.matrix(merged))
}

## older version of yeastParser.  Keep in for a while
##yeastParser <- function(yeastData,
##                        colKeep = c(1, 4, 5, 6, 12, 13)){
##    domain <- read.delim(yeastData, header=FALSE, as.is=TRUE)
##    domain <- domain[,colKeep]
##    colnames(domain) <- c("PROBE", "METHOD", "ID", "IDDES", "INTERPROTMP", "INTERPROTMPDES")
##
##    ##remove version #s from ID
##    domain$"ID" <- sapply(domain$"ID", function(x) unlist(strsplit(x, ".", fixed=TRUE))[1])
##    domain$"INTERPROTMP" <- sapply(domain$"INTERPROTMP", function(x) unlist(strsplit(x, ".", fixed=TRUE))[1])
##
##    ##indeies of methods
##    isDomain <- list()
##    isDomain[["PFAM"]]     <- domain$"METHOD" == "HMMPfam"
##    isDomain[["INTERPRO"]] <- rep(TRUE, length(domain$"METHOD"))
##    
##    merged <- as.matrix(unique(domain$"PROBE"))
##    colnames(merged) <- "PROBE"
##    for(id in names(isDomain)){
##
##        if(id != "INTERPRO"){
##            idindex <- 3
##            desindex <- 4
##        }else{
##            idindex <- 5
##            desindex <- 6
##        }
##
##        tmpDomain <- domain[isDomain[[id]], c(1, idindex, desindex)]
##        ## This one is modified because the PFAM package is ready
##        tmpDomain[,id] <- apply(tmpDomain, 1, function(x) paste(x[3], x[2], sep="@"))
##        
##        temp <- lapply(split(tmpDomain[, id], tmpDomain[, "PROBE"]), unique)
##        temp <- lapply(temp, function(x) x[x!="NULL@NULL"])
##        temp <- sapply(temp, function(x) paste(x, collapse=";"))
##        temp <- as.matrix(temp)
##        temp <- cbind(rownames(temp), temp)
##        rownames(temp) <- NULL
##        colnames(temp) <- c("PROBE", id)
##        merged <- merge(merged, temp, by="PROBE", all.x=TRUE)
##
##        ## doesn't work for some reason
##        ##domain[, "tmp"] <- apply(domain, 1, function(x) ifelse(x[isDomainIndex], paste(x[desindex], x[idindex], sep="@"), "") )
##
##        ##browser()
##    }
##
##    merged[merged==""] <- NA
##    
##    return(as.matrix(merged))
##}
