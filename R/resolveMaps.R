# This function unifies mappings obtained from different source and
# returns a set of unified mappings based on user's specification.
#
# Copyrihgt 2002, J. Zhnag. All rights reserved.
#

resolveMaps <- function(maps, trusted, srcs, colNames = NULL,
                        outName = "", asFile = TRUE){

    doOneRow <- function(vect){
        temp <- vect[!is.null(vect)]
        temp <- temp[!is.na(temp)]
        temp <- temp[toupper(temp) != "NA"]
        if(length(temp) == 0){
            vote <- c(NA, length(vect))
        }
        if(any(!is.na(temp[trusted]))){
            # Get vote from trusted
            vote <- getVote(temp[trusted])
        }else{
            if(any(!is.na(temp[srcs]))){
                # Get vote from sources when there is no vote from trusted
                vote <- getVote(temp[srcs])
            }else{
                # No mapping at all
                vote <- c(NA, length(c(trusted, srcs)))
            }
        }
#        return(vote)

        return(c(vect[setdiff(names(vect), c(trusted, srcs))],
                 vote[1], vote[2]))
    }

    if(!is.null(colNames)){
        colnames(maps) <- colNames
    }
    if(any(!is.element(trusted, colnames(maps)))){
        stop("trusted must be one of the colNames")
    }
    if(asFile && outName == ""){
        outName <- tempfile("tempFile")
    }
#    temp <- apply(maps[c(trusted, srcs)], 1, doOneRow)
#    unified <- cbind(maps[setdiff(names(maps), c(trusted, srcs))], t(temp))
#    if(asFile){
#        write.table(x = unified, file = outName,
#                    quote = FALSE, sep = "\t", row.names = FALSE,
#                    col.names = FALSE)
#        return(outName)
#    }else{
#        return(t(unified))
#    }
    temp <- apply(maps, 1, doOneRow)
    if(asFile){
        write.table(x = t(temp), file = outName,
                    quote = FALSE, sep = "\t", row.names = FALSE,
                    col.names = FALSE)
        return(outName)
    }else{
        return(t(temp))
    }
}

# Finds agreement among a vector containing different sources
getUnified <- function(voters){
    # Find the number of repeatition each element occurs in a vector
    findRep <- function(toFind){
        counts <- length(voters[voters == toFind])
    }
    findSmallest <- function(vec){
        options(warn = -1)
        temp <- try(as.numeric(vec))
        options(warn = 0)
        if(any(is.na(temp))){
            return(sort(vec)[1])
        }else{
            return(sort(temp)[1])
        }
    }

    repeatition <- sapply(voters, findRep)
    majority <- max(repeatition, na.rm = TRUE)
    # If majority = 1, no agreement among sources. Get a mapping by
    # one of the sources
    if(majority == 1){
        #return(c(getNoDup(voters), 1))
      return(c(voters[1], 1))
    }
    # majority <- max(repeatition, na.rm = TRUE)
    found <- c(unique(voters[repeatition == majority]))
    # Only one majority vote
    if(length(found) == 1){
        return(c(found, majority))
    }else{
        return(c(findSmallest(found), majority))
    }
}

# When sources do no agree, get one based on the rules defined here
getNoDup<- function(voters){
    if(length(voters) <= 1){
        return(voters)
    #}else if(length(voters) == 1){
     #   return(voters)
    }else{
        return(voters[1])
        #temp <- voters[1]
        # Find the mapping that has the smallest number of set of characters
        #for(i in 2:length(voters)){
        #    if(temp > voters[i]){
        #        temp <- voters[i]
        #    }
        #}
        #return(temp)
    }
}

hasDelimit <- function(entry, deli = ";"){
    if(any(is.null(entry), is.na(entry))){
        return(FALSE)
    }
    if(gsub(paste(".*(", deli, ").*", sep = ""), "\\1",
                                  as.character(entry)) == ";"){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

getVote <- function(voters, sep = ";"){
    entries <- unlist(sapply(voters, strsplit, sep))
    entries <- entries[!is.null(entries)]
    entries <- entries[!is.na(entries)]
    entries <- entries[toupper(entries) != "NA"]
    # All NAs
    if(length(entries) == 0){
        return(c("NA", length(voters)))
    }
    if(!any(duplicated(entries))){
        # No argeement
        #return(c(getNoDup(entries), 1))
        return(c(entries[1], 1))
    }else{
        # With agreements
        return(getUnified(entries))
    }
}









