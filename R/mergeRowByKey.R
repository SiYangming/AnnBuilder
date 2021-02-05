# This function resolves the one to may relationship beween ids
# contained in a file in one column and the values in another
# column. Duplicating values for the same id will be merged as one
# value separated by ";".
#
# Copyright 2002, J. Zhang. All rights reserved.
#

mergeRowByKey <- function (mergeMe, keyCol = 1, sep = ";"){
    mergeTwoCol <- function(x){
        return(paste(unique(x[, -keyCol]), sep = "", collapse = ";"))
    }
    mergeMultiCol <- function(x){
        return(apply(x, 2, function(y)
                     paste(unique(y), sep = "", collapse = ";")))
    }
    # Returns the original value if mergeMe is a vector
    if(is.null(ncol(mergeMe))){
        return(mergeMe)
    }
    merged <- split.data.frame(mergeMe, factor(mergeMe[, keyCol]))
    if(ncol(mergeMe) == 2){
        merged <- sapply(merged, mergeTwoCol)
        merged <- cbind(names(merged), merged)
        colnames(merged) <- c(colnames(mergeMe)[keyCol],
                              colnames(mergeMe)[-keyCol])
        return(merged)

    }else{
        merged <- sapply(merged, mergeMultiCol)
        merged <- t(merged)
        colnames(merged) <- colnames(mergeMe)
        return(merged)
    }
}




























