# This function writes to an XML file using data contained by a file
# passed either as a file or R object. The file of the R object needs
# to be convertable to a matrix.
#
# Copyright 2002 J. Zhang, all rights reserved.
#

fileToXML <- function (targetName, outName, inName, idColName, colNames,
                       multColNames, typeColNames,  multSep = ";",
                       typeSep = ";", fileSep = "\t", header = FALSE,
                       isFile = TRUE, organism = "human",
                       version = "1.0.0"){

    if(!header && isFile){
        if(is.null(colNames) || is.na(colNames) || colNames == ""){
            stop("Parameter colNames has to be provided!")
        }
    }
    if(isFile){
        fileRead <- as.matrix(read.table(file = inName, header =
                                         header, sep = fileSep,
                                         as.is = TRUE, quote = "",
                                         comment.char = ""))
    }else{
        fileRead <- as.matrix(inName)
    }

    if(colNames != "" && !is.null(colNames) && !is.na(colNames)){
       colnames(fileRead) <- colNames
    }else{
        colNames <- colnames(fileRead)
    }
    # Create the our file and write the header
    file.create(outName)
    writeXMLHeader(outName, colNames, targetName, version, organism)
    # Write the first line of the Data node
    write(x = "<AnnBuilder:Data>", file = outName, append = TRUE)
    # Make starting tags for the Entry nodes
    entries <- paste("<AnnBuilder:Entry id=\"", fileRead[,idColName],
                       "\" describ = \"", idColName, "\">", sep = "")
    # Include Item nodes to Entry nodes using data from all the single
    # value columns
    for(i in setdiff(colNames, c(idColName, multColNames, typeColNames))){
        items <- paste("\t<AnnBuilder:Item name=\"", i,
                       "\" value=\"",
                       .formatStr(as.character(fileRead[,i])),
                       "\" />", sep = "")
        entries <- paste(entries, items, sep = "\n")
    }
    # Process the columns with multiple values separated by multSep
    if(!is.null(multColNames) && !is.null(multColNames) &&
       multColNames != ""){
        for(i in multColNames){
            doMultValue <- function(multCol, sep, name){
                splited <- unlist(strsplit(multCol, sep))
                paste("\t<AnnBuilder:Item name=\"", name, "\" value=\"",
                      .formatStr(as.character(splited)), "\" />", sep = "",
                      collapse = "\n")
            }
            items <- lapply(as.vector(fileRead[,i]), doMultValue, multSep, i)
            entries <- paste(entries, items, sep = "\n")
        }
    }
    # Process the columns with a type attribute
    if(!is.null(typeColNames) && !is.null(typeColNames) &&
       typeColNames != ""){
        for(i in typeColNames){
            doTypeValue <- function(typeCol, sep, name){
                if(typeCol != "" && toupper(typeCol) != "NA" &&
                   !is.na(typeCol) && !is.null(typeCol)){
                    splited <- unlist(lapply(typeCol, twoStepSplit))
                    splited <- cbind(splited, names(splited))
                    paste("\t<AnnBuilder:Item name=\"", name, "\" value=\"",
                          .formatStr(as.character(splited[,1])), "\"",
                          paste(" type=\"", splited[,2], "\"", sep = ""),
                          "/>", sep = "", collapse = "\n")
                }
            }
            items <- lapply(fileRead[,i], doTypeValue, typeSep, i)
            entries <- paste(entries, items, sep = "\n")
        }
    }
    # Collapse entries
    entries <- paste(entries, sep = "",
                     collapse = "\n</AnnBuilder:Entry>\n")
    # Add the rest in
    entries <- paste(entries, "</AnnBuilder:Entry>\n",
                     "</AnnBuilder:Data>\n",
                     "</AnnBuilder:Annotate>\n", sep = "", collapse = "\n")
    #Write to file
    write(entries, file = outName, append = TRUE)

}

.formatStr <- function(tobeDone){
    tobeDone <- gsub("&", "&amp;", tobeDone)
    tobeDone <- gsub("<", "&lt;", tobeDone)
    tobeDone <- gsub(">", "&gt;", tobeDone)
    tobeDone <- gsub("\"", "&quot;", tobeDone)
    tobeDone <- gsub("'", "&apos;", tobeDone)
    return(tobeDone)
}


















