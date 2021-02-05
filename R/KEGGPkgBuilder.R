# This function creates a data package containing mappings between
# KEGG pathway names and pathway ids and enzyme ids and GO ids.
#
# Copyright 2002, J. Zhang. All rights reserved.
#

KEGGPkgBuilder <- function(pkgPath, pkgName = "KEGG", version = "1.0.1",
			  author = list(author = "who",maintainer = "who@email.com")){
    pathwayURL = getKEGGFile("path")
    enzymeURL = getKEGGFile("enzyme")
    geneMapURL = getKEGGFile("geneMap") 
   
    # Create two environments to hold the parsed data
    pNameEnv <- new.env(hash = TRUE, parent = NULL)
    pIdEnv <- new.env(hash = TRUE, parent = NULL)
    makeSrcInfo()
    createEmptyDPkg(pkgName, pkgPath)
    # Write mappings between path ids and path names
    kegg <- KEGG(pathwayURL, fromWeb = TRUE)
    idNName <- getKEGGIDNName(kegg, "")
    #saveMat(as.matrix(cbind(names(, go2All)), pkgName = pkgName,
    #        pkgPath = pkgPath, envName = "ALLLOCUSID", fun = tempFun)
    saveMat(cbind(idNName, names(idNName)), pkgName = pkgName,
            pkgPath = pkgPath, envName = "PATHNAME2ID")
    saveMat(cbind(names(idNName), idNName), pkgName = pkgName,
            pkgPath = pkgPath, envName = "PATHID2NAME")
    # Writes mappings between GO and enzyme
    enzymeIdNGO <- getEIdNName(enzymeURL)
    saveMat(enzymeIdNGO, pkgName = pkgName,
            pkgPath = pkgPath, envName = "ENZYMEID2GO")
    saveMat(getReverseMapping(enzymeIdNGO), pkgName = pkgName,
            pkgPath = pkgPath, envName = "GO2ENZYMEID")
    keggGeneMap <- getKEGGGeneMap(c("Homo sapiens",
                                   "Rattus norvegicus", "Mus musculus",
                                   "Arabidopsis thaliana",
                                    "Drosophila melanogaster",
                                    "Saccharomyces cerevisiae"))
    saveMat(keggGeneMap, pkgName, pkgPath, "EXTID2PATHID", keyCol = 1,
                         valCol = 2, fun = splitEntry)
    saveMat(getReverseMapping(keggGeneMap, sep = ";"), pkgName,
                         pkgPath, "PATHID2EXTID", keyCol = 1,
            valCol = 2, fun = splitEntry)
    repList <- getRepList("KEGG", list(kegg = kegg))
    repList[["PKGNAME"]] <- "KEGG"
    writeDocs("", pkgName, pkgPath, version, author,
               repList, "PKGNAME")
    # Write man pages and other files
    writeDescription(pkgName, pkgPath, version, author,
                      "KEGG", paste("Free for academic use.",
                                   "Non-academic users are requested",
                                   "to\n   obtain a license agreement",
                                   "with KEGG"))
    makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
}


getEIdNName <- function(enzymeURL =
                        "http://www.geneontology.org/external2go/ec2go"){

    kegg <- KEGG(srcUrl = enzymeURL, fromWeb = TRUE)
    fileGot <- readData(kegg)

    getECNGO <- function(vect){
        return(c(vect[1], vect[length(vect)]))
    }
    # Remove the header if any
    if(length(grep("^!", fileGot)) > 0){
        fileGot <- fileGot[-grep("^!", fileGot)]
    }
    # Split entries using a space
    fileGot <- sapply(sapply(fileGot, strsplit, " "), unlist)
    # Create a matrix using the splited data
    fileGot <- matrix(sapply(fileGot, getECNGO), ncol = 2, byrow = TRUE)

    return(fileGot)
}


getKEGGFile <- function(whichOne, organism = "hsa"){
    switch(tolower(whichOne),
           path = return(paste(getSrcUrl("KEGG"), "/",
                    "map_title.tab", sep = "")),
           enzyme = return("http://www.geneontology.org/external2go/ec2go"),
           genemap = return(paste(getSrcUrl("KEGG"),
                       organism, paste(organism, "_gene_map.tab",
                             sep = ""), sep = "/")),
           stop("Invalid option for file name"))
}





