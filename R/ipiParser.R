
ipiParser <- function(ipiData,
                      key="Entrez Gene",
                      acc="AC",
                      target=c("Pfam", "PROSITE"), fromWeb=TRUE) {
    if(fromWeb) {
        tmpFile <- loadFromUrl(ipiData)
        con <- file(tmpFile, "r")
    }else{
        con <- gzfile(ipiData, "rb")
    }
    ipiMatrix <- NULL
    keyDr    <- paste("DR   ", key,    sep="")
    while(TRUE) {
        ans <- readIpiRecord(con)
        keyValue <- getDrValue(ans$record, key)
        accValue <- getDrValue(ans$record, acc)
        if(!is.na(keyValue)) {
            tempVector <- keyValue            
            for(DR in target) {
                drValue <- getDrValue(ans$record, DR)
                drValue <- unlist(strsplit(drValue, ";"))
                drValue <- paste(drValue, accValue, sep="@")
                drValue <- paste(drValue, sep="", collapse=";")
                tempVector <- c(tempVector, drValue)
            }
            
            ipiMatrix <- rbind(ipiMatrix, tempVector)
        }
        if (ans$done) {
            break
        }
    }
    close(con)
    colnames(ipiMatrix) <- c("ENTREZID", toupper(target))
    rownames(ipiMatrix) <- NULL

    merged <- NULL
    for(targetName in target) {
        if(is.null(merged)) {
            merged <- combineDuplicate(ipiMatrix, toupper(targetName))
        } else {
            merged <- merge(merged, combineDuplicate(ipiMatrix, toupper(targetName)), by="ENTREZID", all.x=TRUE)
        }
    }

    return(merged)
}

readIpiRecord <- function(con) {
    buf <- NULL
    atEnd <- FALSE
    done <- FALSE
    while(!done) {
        line <- try(readLines(con, n=1, ok=FALSE), silent=TRUE)
        if (inherits(line, "try-error")) {
            atEnd <- TRUE
            done <- TRUE
        } else if (length(grep("//", line)) > 0) {
            done <- TRUE
        }
        buf <- c(buf, line)
    }
    list(record=buf, done=atEnd)
}

getDrValue <- function(record, DR) {
    switch(DR,
           "Entrez Gene" = findEntrezGene(record),
           "AC" = findAC(record),
           "Pfam" = findPfam(record),
           "PROSITE" = findProsite(record)
           )
}

findEntrezGene <- function(record) {
    value <- getDr(record, "DR   Entrez Gene")
    if(!any(is.na(value))) {
        return(unlist(strsplit(value,   "; ", fixed=TRUE))[2])
    } else {
        return(as.character(NA))
    }
}

findAC <- function(record) {
    value <- record[2]
    result <- unlist(strsplit(value,   "; ", fixed=TRUE))[1]
    result <- sub(";", "", result, fixed=TRUE)
    result <- sub("AC   ", "", result, fixed=TRUE)
    return(result)
}

findPfam <- function(record) {
    value <- getDr(record, "DR   Pfam")
    if(!any(is.na(value))) {
        result <- sapply(strsplit(value,   "; ", fixed=TRUE), function(x) x[2])
        return(paste(result, sep = "", collapse = ";"))
    } else {
        return(as.character(NA))
    }
}

findProsite <- function(record) {
    value <- getDr(record, "DR   PROSITE")
    if(!any(is.na(value))) {
        result <- sapply(strsplit(value,   "; ", fixed=TRUE), function(x) x[2])
        return(paste(result, sep = "", collapse = ";"))
    } else {
        return(as.character(NA))
    }
}

parseIpiRecord <- function(rec) {
    print("parsing record: ", substr(rec, 1, 10))
}

getDr <- function(record, DR) {
    indexDR <- grep(DR, record, value=TRUE)
    if(length(indexDR)>0) {
        return(indexDR)
    } else {
        return(as.character(NA))
    }
}

combineDuplicate <- function(ipiMatrix, targetName, byName="ENTREZID") {
    ipiSplit <- split(ipiMatrix[,targetName], ipiMatrix[,byName])
    ipiSplit <- sapply(ipiSplit, function(x) paste(x, collapse=";"))
    ipiSplit <- as.matrix(ipiSplit)
    ipiSplit <- cbind(rownames(ipiSplit), ipiSplit)
    rownames(ipiSplit) <- NULL
    colnames(ipiSplit) <- c(byName, targetName)
    return(ipiSplit)
}

