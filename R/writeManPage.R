# This function writes all the elements other than rda files required
# for a data package
writeAccessory <- function(pkgName, pkgPath, organism, version,
                         author = list(author = "who",
                         maintainer = "My Name <who@email.net>"),
                         dataSrc = "public data repositories",
                         license = "The Artistic License, Version 2.0"){
    # Write to DESCRIPTION
    writeDescription(pkgName, pkgPath, version, author, dataSrc, license)
    # Write Rd files for all the rda files that already exist in data
    AnnInfo <- as.list(get("AnnInfo", env = .GlobalEnv))
    for(i in getAllRdaName(pkgName, pkgPath)){
        writeManPage(pkgName, pkgPath, i, organism,
                     src = get(i, AnnInfo)$src)
    }
    # Write a function that displays data information when called
    # using pkgName() and the man page for the function
    writeFun (pkgPath, pkgName, organism = organism)
    writeMan4Fun(pkgName, pkgPath, organism)
    # Write .First.lib
    writeZZZ(pkgPath, pkgName)
    # Put the pretty print function for QC data and man page in
#    putPrettyPrint(pkgName, pkgPath)
}

# Functions that write the documentation files for a data package.
#
# Copyright 2002, J. Zhang, all rights reserved.
#

writeManPage <- function(pkgName, pkgPath, manName, organism = "human",
                         src = "ll", isEnv = TRUE){

    path <- file.path(pkgPath, pkgName, "man")
    filename <- file.path(path, paste(pkgName, formatName(manName),
                                    ".Rd", sep = ""))
    if(is.na(src)){
        srcBuild <- "N/A"
    }else{
        if(toupper(organism) == "YEAST"){
            srcBuild <- getBuild4Yeast(src, manName)
        }else{
            srcBuild <- getSrcBuiltNRef(src, organism)
        }
    }
    AnnInfo <- get("AnnInfo")
    write(paste("\\name{", paste(pkgName, manName, sep = ""), "}\n",
                "\\alias{", paste(pkgName, manName, sep = ""), "}\n",
                "\\title{", "An annotation data file for ", manName,
                " in the ", pkgName, " package}\n",
                "\\description{\n", get(manName, AnnInfo)$long, "\n}\n",
                srcBuild, getExample(pkgName, manName, isEnv),
                "\\keyword{datasets}\n", sep = ""), filename)
}

getExample <- function(pkgName, manName, isEnv = TRUE){
    if(isEnv){
        return(paste("\\examples{\n", paste("\trequire(\"annotate\")",
                                            "|| stop(\"annotate",
                                            "unavailable\")\n"),
                     paste(paste("\txx <- ls(env = ",
                           paste(pkgName, formatName(manName), sep = ""),
                           ")", sep = ""),
                           "\tif(length(xx) > 0){",
                           "\t\t# Using get for value of the first key",
                           paste("\t\tget(xx[1],", paste(pkgName,
                              formatName(manName), sep = ""), ")"),
                           "\t\t#Using mget for a few keys",
                           "\t\tif(length(xx) >= 3){",
                           paste("\t\t\tmget(xx[1:3],", paste(pkgName,
                              formatName(manName), sep = ""), ")"),
                           "\t\t\t#Using lookUp of annotate(> 1.3.4)",
                           paste("\t\t\tlookUp(xx[1:3],\"", pkgName, "\",\"",
                                 manName, "\")", sep = ""), "\t\t}",
                           sep = "\n"),
                           "\n\t}\n}\n", sep = ""))
    }else{
        return(paste("\\examples{\n",
               paste("length(", pkgName, manName, ")", sep = ""),
                     "\n}\n", sep = ""))
    }
}

getSrcBuiltNRef <- function(src, organism){
    switch(toupper(src),
           "LL" = dataSource <- "LocusLink",
           "UG" = dataSource <- "UniGene",
           "KEGG" = dataSource <- "KEGG",
           "GO" = dataSource <- "Gene Ontology",
           "GP" = dataSource <- "Human Genome Project",
           "YG" = dataSource <- "Yeast Genome Project",
           "HG" = dataSource <- "HomoloGene",
           "REFSEQ" = dataSource <- "Reference Sequence Project",
           dataSource <- src)

    if(src == ""){
        return(paste("\\details{\n",
               paste(paste("Package built:", date())), "\n}\n", sep = ""))
    }else{
        return(paste("\\details{\n", paste("Mappings were based on data ",
                                           "provided by ", dataSource,
                                           "\n\n", sep = ""),
                     paste(paste("Source data built:",
                                 getSrcNBuilt(src, organism), "\n",
                                 paste("Package built:", date()),
                                 sep = "")), "\n}\n",
                     "\\references{\n", paste("\\url{", gsub("_", "\\\\_", getSrcUrl(src)), "}",
                                              sep = ""), "\n}\n", sep = ""))
    }
}

# This function writes the man page for the function that displays
# data information when invoked by typing pkgName()
writeMan4Fun <- function(pkgName, pkgPath, organism = "human", QCList,
                         dSrc = "all"){
    fileName <- file.path(pkgPath, pkgName, "man",
                          paste(pkgName, ".Rd", sep = ""))
    if(dSrc == "all"){
        dSrc <- getDSrc(organism)
    }
    item <- getSrcNBuilt(dSrc, organism)
    item <- paste(item, sep = "\n")
    write(paste("\\name{", pkgName, "}\n",
                "\\alias{", pkgName, "}\n",
                "\\title{",
                "Genomic Annotation data package built with AnnBuilder}\n",
                "\\description{\nThe package is built using a ",
                "downloadable R package - AnnBuilder (download and build ",
                "your own) from www.bioconductor.org using the ",
                "following public data sources:\n", item,
                "\n\nThe function ", pkgName, "() provides ",
                "information about the binary data files}\n",
                "\\keyword{datasets}\n", sep = ""), fileName)

}

escapeLatexChr <- function(item){
    item <- gsub("_", "\\\\_", item)
    item <- gsub("#", "\\\\#", item)
    item <- gsub("&", "\\\\&", item)
    item <- gsub("%", "\\\\%", item)
    item <- gsub("\\{", "\\\\{", item)
    item <- gsub("\\}", "\\\\}", item)
    return(item)
}

getSrcNBuilt <- function(dSrc, organism){
    items <- ""
    for(i in dSrc){
        switch(toupper(i),
               "LL" = toAdd <- paste("LocusLink",
                                      getUrlNBuilt("LL", organism)),
               "GP" = toAdd <- paste("Human Genome Project",
                                 getUrlNBuilt("GP", organism)),
               "UG" = toAdd <- paste("UniGene",
                                      getUrlNBuilt("UG", organism)),
               "GO" = toAdd <- paste("Gene Ontology Consortium",
                                      getUrlNBuilt("GO", organism)),
               "KEGG" = toAdd <- paste("KEGG",
                                      getUrlNBuilt("KEGG", organism)),
               "YG" = toAdd <- paste("Yeast Genome",
                                      getUrlNBuilt("YG", organism)),
               "HG" = toAdd <- paste("Homology",
                                      getUrlNBuilt("HG"), organism),
               "REFSEQ" = toAdd <- paste("Reference Sequence project",
                                      getUrlNBuilt("REFSEQ", organism),
                                      organism),
               warning(paste("Source", i, "is not a valid name!")))
        items <- paste(items, "\n", toAdd, sep = "")
    }
    return(items)
}

getUrlNBuilt <- function(src, organism){
    paste("built: ", escapeLatexChr(getSrcBuilt(src, organism)), ".",
          "\\url{", escapeLatexChr(getSrcUrl(src, organism)), "}.", sep = "")
}

#getItem <- function(name, content){
#    paste("\\item{", name, "}{",content, "}", sep = "")
#}

formatName <- function (toFormat){
    gsub("\\.*(_)*", "", toFormat)
}

writeREADME <- function(pkgPath, pkgName, urls){

    OPENING <- paste("The binary data files in this package were constructed",
                     "using data from the following sources identified",
                     "by their URLs:")
    CLOSING <- paste("We highly appreciate their efforts of making",
                     "the data available publically")


    fileName <- file.path(pkgPath, pkgName, "README")

    write(OPENING, file = fileName, append = FALSE)
    for(i in urls)
        write(i, file = fileName, append = TRUE)
    write(CLOSING, file = fileName, append = TRUE)
}


writeDescription <- function(pkgName, pkgPath, version, author,
                             dataSrc = "public data repositories",
                             license = "The Artistic License, Version 2.0") {
    path <- file.path(pkgPath, pkgName)
    fileName <- file.path(path, "DESCRIPTION")
    if(file.exists(fileName)){
        unlink(fileName)
    }
    rVersion <- paste(R.Version()[c("major", "minor")], collapse=".")
    ## read destriptionInfo.txt
    data("descriptionInfo")
    descriptionInfo <- get("descriptionInfo")
    slots <- descriptionInfo[which(descriptionInfo[,"biocPkgName"]==pkgName),]
    
    write(paste("Package:", pkgName), file=fileName, append=FALSE)
    ## nrow(slots)==0 means that there is no such pkgName in the descriptionInfo
    if(nrow(slots)==0 || slots["chipName"]=="") {
        write(paste("Title: A data package containing",
                    "annotation data for", pkgName),
              file=fileName, append=TRUE)
        write(paste("Description: Annotation data file for",
                    pkgName, "assembled using data\n   from",
                    dataSrc),
              file=fileName, append=TRUE)
    }else{
        write(paste("Title: ", slots["manufacturer"], " ", slots["chipName"], " Annotation Data (", pkgName, ")", sep=""),
              file=fileName, append=TRUE)
        write(paste("Description: ", slots["manufacturer"], " ", slots["chipName"], " annotation data (", pkgName, 
                    ") assembled using data from ", dataSrc, sep=""),
              file=fileName, append=TRUE)
    }
    write(paste("Version:", version,
                "\nCreated:", date(),
                "\nAuthor:", paste(author[["authors"]], collapse=", "),
                "\nMaintainer:", author[["maintainer"]],
                "\nLazyData: yes",
                paste("\nDepends: R(>= ", rVersion, ")", sep=""),
                "\nLicense:", license),
          file=fileName, append=TRUE)
    if(pkgName=="GO") {
        write("Suggests: annotate",
          file=fileName, append=TRUE)
    }
    if(nrow(slots)!=0) {
       	slots <- slots[,-which(colnames(slots)=="biocPkgName")]
       	slots[is.na(slots)] <- ""
	slotNames <- names(slots)
	if (slots["chipName"]=="") {## YEAST, *LLMappings, etc
        	for(i in 1:length(slots)) {
			if (slots[i]!="") {
            			write(paste(slotNames[i], ": ", slots[i], sep=""),
        	      		file=fileName, append=TRUE)
			}
        	}
	} else {## chip annotation packages
        	for(i in 1:length(slots)) {
            		write(paste(slotNames[i], ": ", slots[i], sep=""),
        	      	file=fileName, append=TRUE)
        	}
	}
    }  
}

getDSrc <- function(organism){
    switch(tolower(organism),
           "human" = ,
           "rat" = ,
           "mouse" = return(c("LL", "GO", "KEGG", "GP", "UG")),
           "yeast" = return(c("KEGG", "GO", "YG")),
           "kegg" = return("KEGG"),
           "go" = return("GO"),
           "hg" = return("HG"),
           stop(paste("Organism", organism, "not supported!")))
}

# This function searches the package to be built for names of the rda
# files already created
getAllRdaName <- function(pkgName, pkgPath){
    rdaNames <- gsub(".rda", "",
                     list.files(file.path(pkgPath, pkgName, "data")))
    rdaNames <- sub(pkgName, "", rdaNames)
    return(rdaNames[!is.element(rdaNames, c("QC", "print.ABQCList"))])
}

writeWarning <- function(pkgName, oldObj, newObj, useMsg=F) {
    if (useMsg) {
    		toWrite <- paste(
                     pkgName, oldObj, "func <- function() {\n",
                     "\t.Deprecated(msg=\"", newObj, "\")\n",
                     "\t", pkgName, oldObj, "_DEPRECATED\n",
                     "}\n",
                     "makeActiveBinding(\"", pkgName, oldObj, "\", ",
                     pkgName, oldObj, "func, ns)\n",
                     sep="", collapse="")
    } else { 
    		toWrite <- paste(
                     pkgName, oldObj, "func <- function() {\n",
                     "\t", pkgName, oldObj, " <- function(){\n",
                     "\t\t.Deprecated(\"", pkgName, newObj, "\")\n",
                     "\t}\n",
                     "\t", pkgName, oldObj, "()\n",
                     "\t", pkgName, newObj, "\n",
                     "}\n",
                     "makeActiveBinding(\"", pkgName, oldObj, "\", ",
                     pkgName, oldObj, "func, ns)\n",
                     sep="", collapse="")
   }
   toWrite
}

writeFun <- function (pkgPath, pkgName, organism = "human"){
    Rdir <- file.path(pkgPath, pkgName, "R")
    fileName <- file.path(Rdir, paste(pkgName, ".R", sep = ""))
    qc <- paste(pkgName, "QC", sep="")
    toWrite <- paste(pkgName, " <- function() cat(", qc, ")\n",
                     ".no_S3_generics = TRUE\n",
                     "ns <- asNamespace(\"", pkgName, "\")\n",
                     sep="", collapse="") 
    specialPkgs <- c("cMAP", "KEGG", "humanLLMappings", "mouseLLMappings", "ratLLMappings", "YEAST", "yeast2", "ygs98", "PFAM")
    if (length(grep(".*CHRLOC$", pkgName))>0) {
	toWrite <- paste(toWrite,
                     writeWarning(pkgName, "LOCUSID2CHR", "ENTREZID2CHR"),
                     sep="", collapse="")
	depObjs <- paste(pkgName, "LOCUSID2CHR", sep="", collapse="")
    } else if (pkgName=="GO"){
	toWrite <- paste(toWrite,
                     writeWarning(pkgName, "LOCUSID", "ENTREZID"),
                     writeWarning(pkgName, "LOCUSID2GO", "ENTREZID2GO"),
                     writeWarning(pkgName, "ALLLOCUSID", "ALLENTREZID"),
                     sep="")
	depObjs <- paste(pkgName, "LOCUSID, ", pkgName, "LOCUSID2GO, ", 
                     pkgName, "ALLLOCUSID", sep="", collapse="")
    } else if (pkgName=="ath1121501"||pkgName=="ag"){
        toWrite <- paste(toWrite,
                     "delayedAssign(\"", pkgName, "ENTREZID\", ",
                     pkgName, "ACCNUM)\n",
                     writeWarning(pkgName, "LOCUSID", "ENTREZID"),
                     sep="")
	depObjs <- paste(pkgName, "LOCUSID, ", pkgName, "ENTREZID", 
                     sep="", collapse="")
    } else if (! pkgName %in% specialPkgs){
## 	toWrite <- paste(toWrite,
##                      writeWarning(pkgName, "LOCUSID", "ENTREZID"),
##                      #writeWarning(pkgName, "SUMFUNC", "This environment never contains any data, and will be removed in the next release.", useMsg=T),
##                      sep="", collapse="") 
## 	depObjs <- paste(pkgName, "LOCUSID", sep="", collapse=", ")
## 	#depObjs <- paste(pkgName, c("LOCUSID", "SUMFUNC"), sep="", collapse=", ")
    }
    write(toWrite, file = fileName)
    nsFile <- file.path(pkgPath, pkgName, "NAMESPACE")
    if (pkgName %in% specialPkgs) 
	toWrite <- paste("export(", pkgName, ")\n", sep="", collapse="")
    else
##     	toWrite <- paste("export(", pkgName, ", ", depObjs, ")\n", sep="", collapse="")
    write(toWrite, file=nsFile)
}

writeZZZ <- function(pkgPath, pkgName){
    return(TRUE)
}


writeMan4QC <- function(pkgName, pkgPath){

    path <- file.path(pkgPath, pkgName, "man")
    filename <- file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "QC.Rd", sep = ""))
    write(paste("\\name{", pkgName, "QC}\n",
                "\\alias{", pkgName, "QC}\n",
                "\\title{", "Quality control information for ",
                pkgName, "}\n",
                "\\description{\n A data file containing statistics ",
                "for all the data files in this package. The data ",
                "can be used for quality control purpose\n}\n",
                "\\details{\n", paste("This file contains quality control ",
                            "information that can be displayed by ",
                            "typing ", pkgName, "() after loading ",
                            "the package using library(", pkgName,
                            "). The follow items are included:\n\n",
                            "Date built - Date when the package was built.",
                            "\n\nNumber of probes - total number of ",
                            "probes included\n\n", "Probe number ",
                            "missmatch - if the total number of probes",
                            " of any of the data file is different ",
                            "from a base file used to check the data files",
                            " the name of the data file will be listed\n\n",
                            "Probe missmatch - if any of probes in a data",
                            " file missmatched that of the base file, ",
                            " the name of the data file will be listed\n\n",
                            " Mappings found for probe based files - ",
                            " number of mappings obtained for the total ",
                            "number of probes\n\n", "Mappings found for ",
                            "non-probe based files - total number of ",
                            "mappings obtained", sep = ""), "\n}\n",
                "\\keyword{datasets}\n", sep = ""), filename)
}

# Get build information for yeast genome source data
getBuild4Yeast <- function(src, manName){
    if(is.element(manName, c("GENENAME", "PMID", "GO", "ORF", "CHR",
                             "CHRLOC", "DESCRIPTION", "CHRLENGTHS",
                             "GO2PROBE", "GO2ALLPROBES", "ALIAS",
                             "#CHROMSEQ#"))){
        if(!any(manName == c("DESCRIPTION", "#CHROMSEQ#"))){
            probeDesc <- " PROBE ids are Open Reading Frame (ORF) ids."
        }else{
            probeDesc <- ""
        }
        return(paste("\\details{\n", "Annotation based on data ",
                     "provided by Yeast Genome project.", probeDesc, "\n\n",
                     paste(paste("Source data built:",
                     "Yeast Genome data are built at various time ",
                     "intervals. Sources used were downloaded ",
                     date(), "\n",
                     paste("Package built:", date()), sep = "")), "\n}\n",
                     "\\references{\n", paste("\\url{", getYGUrl(), "}",
                                              sep = ""), "\n}\n", sep = ""))
    }else{
        return(getSrcBuiltNRef(src, "YG"))
    }
}

