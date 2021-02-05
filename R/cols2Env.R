# This function takes two data columns from a data file (e. g. two
# columns) and generates an environment object with
# values in one of the columns as keys and values in the other as
# values. A given row in either column may have multiple values
# separated by a separator(e. g. "a;b").
#
# Copyright 2002, J. Zhang, all rights reserved
#

cols2Env <- function(cols, colNames, keyColName = colNames[1], sep = ";"){

    # Environment to return
    dataEnv <- new.env(hash = TRUE, parent = emptyenv())

    # Force a matrix
    cols <- as.matrix(cols)

    if(ncol(cols) != 2){
        stop(paste("The input data", cols, "do not have two columns"))
    }
    colnames(cols) <- colNames
    # Determine the column name
    valColName <- colNames[colNames != keyColName]
    # Match cols by finding matches for all the potential multiple
    # matches
    cols <- matchAll(cols, keyColName)
    if(is.null(cols)){
        return(dataEnv)
    }
    colnames(cols) <- colNames
    #Get unduplicated row
    unDups <- cols[!duplicated(cols[,keyColName]),]
    if(length(unDups) != 0){
        # In case there is only one entry
        if(is.null(nrow(unDups))){
            unDups <- matrix(unDups, ncol = 2, byrow = TRUE)
        }
        colnames(unDups) <- colNames
        # Write unique keys and values
        multiassign(unDups[, keyColName],
                lapply(unDups[, valColName], splitEntry, sep = sep), dataEnv)
    }
    # Process the duplications
    dups <- cols[duplicated(cols[,keyColName]),]
    if(length(dups) != 0){
        if(is.null(nrow(dups))){
            dups <- matrix(dups, ncol = 2, byrow = TRUE)
        }
        colnames(dups) <- colNames
        for(i in 1:nrow(dups)){
            assign(dups[i,keyColName], c(get(dups[i,
                       keyColName], dataEnv), unique(dups[i,valColName])),
                   dataEnv)
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
