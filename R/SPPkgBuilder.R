SPPkgBuilder <- function(pkgPath, version, author, fromWeb = TRUE,
  url = "ftp://ftp.ebi.ac.uk/pub/databases/swissprot/release/sprot41.dat"){

    dataSrc <- paste("Swiss-Prot (\\url{", url, "})", sep = "")
    # Just use LL to create a Swiss-Prot object
    sp <- LL(srcUrl = url, parser = file.path(.path.package("AnnBuilder"),
                           "data", "SPParser"))
    spData <- as.matrix(parseData(sp, ncol = 25, fromWeb = fromWeb))
    colnames(spData) <- getEnvNames()
    spData <- spData[2:nrow(spData),]
    spData[spData == "NA"] <- NA
    # Write AnnInfo to .GlobalEnv to make Annotation information available
    makeSrcInfo()
    AnnInfo <- as.list(get("AnnInfo", env = .GlobalEnv))
    for(i in c("human", "mouse", "rat")){
        pkgName <- paste(i, "SP", sep = "")
        tempData <- spData[spData[, "ORGANISM"] == toupper(i),]
        # Create a data package with no data
        createEmptyDPkg(pkgName, pkgPath, force = TRUE)
        # Write the template for man pages
        writeManPage(pkgName, pkgPath, "#ONETOONEMAN#", "",
                     src = get("#ONETOONEMAN#", AnnInfo)$src)
        writeManPage(pkgName, pkgPath, "#ONETOMANYMAN#", "",
                     src = get("#ONETOMANYMAN#", AnnInfo)$src)
        for(j in colnames(tempData)[!is.element(colnames(tempData),
                                                c("AC", "ORGANISM"))]){
            if(isOneToOne(j)){
                # AC may have multiple values. Collapse before writing
                saveData2Env(getReverseMapping(cbind(tempData[, j],
                                   gsub(" +", "", tempData[, "AC"]))),
                             fun = function(x) x, pkgName, pkgPath, j)
                copySubstitute(file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "#ONETOONEMAN#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, j, ".Rd", sep = "")),
                   list(ONETOONEMAN = j,
                        VALUE = tolower(getDetailV(j)),
                        KEY = "Swiss-Prot accession number",
                        ENVNAME = paste(pkgName, j, sep = ""),
                        DATASOURCE = dataSrc, ), "#")
            }else{
                saveData2Env(getReverseMapping(cbind(
                                   gsub(" +", "", tempData[, j]),
                             gsub(" +", "", tempData[, "AC"]))),
                             twoStepSplit,
                             pkgName, pkgPath, j)
                copySubstitute(file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "#ONETOMANYMAN#.Rd", sep = "")),
                           file.path(pkgPath, pkgName,
                           "man", paste(pkgName, j, ".Rd", sep = "")),
                   list(ONETOMANYMAN = j,
                        VALUE = tolower(getDetailV(j)),
                        KEY = "protein accession number",
                        ENVNAME = paste(pkgName, j, sep = ""),
                        DATASOURCE = dataSrc, ), "#")
            }
        }
        # Do reverse mapping between MIM ids and accession number
        AC2MIM <- tempData[!is.na(tempData[, "MIM"]), c("AC", "MIM")]
        AC2MIM <- cbind(gsub(" +", "", AC2MIM[,1]), AC2MIM[,2])
        saveData2Env(getReverseMapping(AC2MIM),
                              twoStepSplit, pkgName, pkgPath, "MIM2AC")
        copySubstitute(file.path(pkgPath, pkgName, "man",
                       paste(pkgName, "#ONETOMANYMAN#.Rd", sep = "")),
                       file.path(pkgPath, pkgName, "man",
                              paste(pkgName, "MIM2AC", ".Rd", sep = "")),
                       list(ONETOMANYMAN = "MIM2AC",
                            VALUE = "Swess-Prot protein accession number",
                            KEY = "MIM id",
                            ENVNAME = paste(pkgName, "MIM2AC", sep = ""),
                            DATASOURCE = dataSrc, ), "#")
        # Write man pages and so on
        writeDescription(pkgName, pkgPath, version, author,
                         dataSrc = "NCBI")
        writeFun (pkgPath, pkgName, organism = organism)
        writeMan4Fun(pkgName, pkgPath)
        # Write .First.lib
        writeZZZ(pkgPath, pkgName)
        # Write the quality control data
        getDPStats("", pkgName, pkgPath)
        writeMan4QC(pkgName, pkgPath)
        file.remove(file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "#ONETOONEMAN#.Rd", sep = "")))
        file.remove(file.path(pkgPath, pkgName, "man",
                           paste(pkgName, "#ONETOMANYMAN#.Rd", sep = "")))
        makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
    }
}

getDetailV <- function(key){
    values <- c("Swiss-Prot accession number","protein name",
                "specise name", "type of molecule",
                "total number of aminoacids", "PubMed id",
                paste("description of related protein sequence(s)",
                      "produced by alternative splicing or initiation",
                      "codons"),
                "description of reactions catalyzed",
                "possible errors and/or grounds for confusion",
                "an enzyme cofactor",
                paste("description of developmentally specific",
                      "expression of a protein"),
                paste("description of the disease(s) associated with",
                      "a deficiency of a protein"),
                "description of the domain structure of a protein",
                "description of an enzyme regulatory mechanism",
                "description of the function(s) of a protein",
                paste("description of the compound(s) which stimulate",
                      "the synthesis of a protein"),
                paste("weight of a protein or part of a protein as",
                      "determined by mass spectrometric methods"),
                paste("descriptionof the use of a protein as a",
                      "pharmaceutical drug"),
                "description of polymorphism(s)",
                "description of a posttranslational modification",
                paste("description of the similarities(s) of a",
                      "protein with orther proteins"),
                paste("description of the subcellular lacation of",
                      "the mature protein"),
                "desction of the quaternary structure of a protein",
                "description of the tissue specificity of a protein",
                "MIM id")
    names(values) <- getEnvNames()
    return(values[key])
}

getEnvNames <- function(){
    return(toupper(c("ac", "name", "organism", "type", "length", "pmid",
               "AP", "CA", "caution", "cofactor", "DS", "disease",
               "domain", "ER", "function", "induction", "MS", "pharm",
               "polym", "PTM", "similarity", "SL", "subunit", "TS","MIM")))
}

isOneToOne <- function(envName){
    if(any(envName == c("PMID"))){
        return(FALSE)
    }else{
        return(TRUE)
    }
}

