cMapPathBuilder <- function(cartaName, keggName, pkgName = "cMAP",
                            pkgPath,  version = "1.1.0",
                            author = list(author = "anonymous",
                            maintainer = "anonymous@email.com"), lazyLoad = TRUE){
    biocarta <- cMAPParser(cartaName)
    kegg <- cMAPParser(keggName)

    createEmptyDPkg(pkgName = pkgName, pkgPath, force = TRUE)
    for(i in names(biocarta)){
        saveList(biocarta[[i]], pkgName, pkgPath,
                                 paste("CARTA", toupper(i), sep = ""))
    }
    for(i in names(kegg)){
        saveList(kegg[[i]], pkgName, pkgPath,
                                 paste("KEGG", toupper(i), sep = ""))
    }
    repList <- list(CMAPSOURCE = paste("cMAP: \\\\url{http://cmap.nci.",
                    "nih.gov/PW/Download}. Build: Unavailable. Downloaded:",
                    date(), sep = ""), DATE = date())
    repList[["PKGNAME"]] <- pkgName

    for(i in c("CARTA", "KEGG")){
        repList[["DSOURCE"]] <- i
        for(j in c("MOLECULE", "INTERACTION", "PATHWAY")){
            rdName <- file.path(pkgPath, pkgName, "man",
                            paste(pkgName, i, j, ".Rd", sep = ""))
            copied <- try(copySubstitute(file.path(
                          .path.package("AnnBuilder"), "templates",
                          paste("PKGNAME", j, ".Rd",  sep = "")),
                                     rdName, repList, "#"))
        }
    }
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates",
                             "GENERAL.Rd"), file.path(pkgPath, pkgName,
                             "man", paste(pkgName, ".Rd", sep = "")),
                              list(PKGNAME = pkgName, SOURCENBUILT =
                                   paste("\n\t", repList, sep = "",
                                         collapse = "\n")), "#")
    writeDescription(pkgName, pkgPath, version, author)
    writeZZZ(pkgPath, pkgName)
    writeFun(pkgPath, pkgName)
    getDPStats("", pkgName, pkgPath)
    writeMan4QC(pkgName, pkgPath)
    copySubstitute(file.path(.path.package("AnnBuilder"),
                        "templates", "PKGNAMEQC.Rd"),
                        file.path(pkgPath, pkgName, "man",
                        paste(pkgName, "QC.Rd", sep = "")),
                        list(PKGNAME = pkgName), "#")
    writeDatalist(pkgName, pkgPath)
    if(lazyLoad){
      makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
    }
}

# Parses cMAP pathway interaction data
# sourceFile should have been downloaded from
# http://cmap.nci.nih.gov/PW/Download

cMAPParser <- function(sourceFile) {
    # Keeps a temporary value
    isSource <- FALSE
    isORG <- FALSE
    isName <- FALSE
    isSName <- FALSE
    isInter <- FALSE
    isPath <- FALSE
    isMol <- TRUE

    molList <- list()
    mID <- NA
    labelType <- NA
    nameType <- NA
    molType <- NA
    compId <- NA
    component <- list()
    member <- NA
    firstID <- NA

    interList <- list()
    interID <- NA
    process <- NA
    interComp <- list()
    edge <- NA
    func <- NA
    location <- NA
    activity <- NA
    interCompId <- NA
    condition <- NA
    reverse <- NA
    role <- NA
    dSource <- NA

    pathList <- list()
    pathId <- NA
    organism <- NA
    pathName <- NA
    shortName <- NA
    pathComp <- NA

    writeInter <- function(){
        #print("Before")
        interList[[strsplit(interID,"@")[[1]][1]]] <<- list(source = dSource, process = strsplit(interID,"@")[[1]][2],
                                     reversible = reverse,
                                     condition = condition,
                                     component = interComp)
        print(paste("writeIntercomponetn",condition ,sep=" "))
        interID <<- NA
        process <<- NA
        interComp <<- list()
        edge <<- NA
        func <<- NA
        location <<- NA
        activity <<- NA
        interCompId <<- NA
        condition <<- NA
        reverse <<- NA
        dSource <<- NA
        #print("after")
    }

    writeInterComp <- function(){
        interComp[[length(interComp) + 1]] <<- list(id = strsplit(interCompId,"@")[[1]][1],
                                              edge = strsplit(interCompId,"@")[[1]][2], 
                                              location = location,
                                              activity = activity)
        interCompId <<- NA
        edge <<- NA
        location <<- NA
        activity <<- NA
    }
    writeMol <- function(){
        #print(paste("before", mID))
        mID <- strsplit(molType,"@")[[1]][1]
        molList[[mID]] <<- list(type =  strsplit(molType,"@")[[1]][2],
                                extid = twoStepSplit(paste(unique(nameType),
                                                      collapse = ";")),
                                component = component,
                                member = twoStepSplit(paste(unique(member),
                                      collapse = ";"), asNumeric = TRUE))
        #print("after")
        reset()
    }
    doLabelType <- function(attr){
        #print("before")
        switch(attr["label_type"],
               "location" = location <<- as.character(attr["value"]),
               #"process-type" = process <<- as.character(attr["value"]),
               "reversible" = reverse <<- as.logical(
                              ifelse(attr["value"] == "no", FALSE, TRUE)),
               #"process-condition" = doCondition(attr),
               #"edge-type" = edge <<- as.character(attr["value"]),
               "activity-state" = activity <<- as.character(attr["value"]),
               warning("Unknown attribute name"))
        #print(paste("doLabelType",location,sep=" "))
    }
    doCondition <- function(attr){
        if(is.na(condition[1])){
            condition <<- attr["value"]
        }else{
            condition <<- c(condition, attr["value"])
        }
    }
    doLocation <- function(attr){
        if(isMol){
            if(length(component) == 1){
                names(component) <<- as.character(attr["value"])
            }else{
                names(component) <<- c(names(component)[-length(component)],
                                       as.character(attr["value"]))
            }
        }else{
            location <<- as.character(attr["value"])
        }
    }
    writeComp <- function(){
        component[[length(component) + 1]] <<- list(id = compId,
                                              location = location,
                                              activity = activity)
        compId <<- NA
        location <<- NA
        activity <<- NA
    }
    doNameType <- function(attr){
        #print(paste("before", attr["value"]))
        if(!isPath){
            if(!is.na(nameType[1])){
                nameType <<- c(nameType, paste(attr["value"],
                                           attr["name_type"], sep = "@"))
            }else{
                nameType <<-  paste(attr["value"],
                                    attr["name_type"], sep = "@")
            }
        }else{
            isName <<- TRUE
        }
        #print("after")
    }
    doMolType <- function(attr){
        #print("before")
            if(!is.na(molType[1])){
                molType <<- c(molType, paste(attr["id"], attr["molecule_type"], sep = "@"))
            }
	    else{
                 molType <<- paste(attr["id"], attr["molecule_type"], sep = "@")
            }
        #print("after")
    }

    doComplex <- function(attr){
        if(is.na(component[1])){
            component <<- as.character(attr["idref"])
            names(component) <<- NA
        }else{
            component <<- c(component, attr["idref"])
            names(component) <<- c(names(component)[-length(component)],
                                   NA)
        }
    }
    
    doFamily <- function(attr){
        #print("before")
        if(is.na(member[1])){
            member <<- paste(attr["member_molecule_id"],
                             attr["family_molecule_id"], sep = "@")
        }else{
            member <<- c(member,  paste(attr["member_molecule_id"],
                             attr["family_molecule_id"], sep = "@"))
        }
        #print("after")
    }
    reset <- function(){
         mID <<- NA
         labelType <<- NA
         nameType <<- NA
         molType <<- NA
         component <<- list()
         member <<- NA
    }

    writePath <- function(){
        pathList[[shortName]] <<- list(id = pathId, organism = organism, source = dSource,
                              name = pathName, component = unique(pathComp))
        
        print (paste(dSource, pathId, organism, pathName, sep=" "))                 
        dSource <<- NA
        pathId <<- NA
        organism <<- NA
        pathName <<- NA
        shortName <<- NA
        pathComp <<- NA
    }

    doPathComp <- function(attr){
        #print("before")
        if(is.na(pathComp[1])){
            pathComp <<- as.integer(attr["idref"])
        }else{
            pathComp <<- c(pathComp, as.integer(attr["idref"]))
        }
        #print("after")
    }

    text <- function(x){
        if(isSource){
            dSource <<- x
            isSource <<- FALSE
        }else if(isName){
            pathName <<- x
            isName <<- FALSE
        }else if(isSName){
            shortName <<- x
            isSName <<- FALSE
        }else if(isORG){
            #print(paste("organism = ", x))
            organism <<- x
            isORG <<- FALSE
        }
    }

    # Event handlers for the parser
    getHandler <- function(){
        # Gets the value for desired element
        startElement <- function(name, attrs){
           # print(paste(name, ";", attrs))
            switch(name,
                   "Molecule" = doMolType(attrs),
                   "Label" = doLabelType(attrs),
                   "Name" = doNameType(attrs),
                   "Condition" = condition <<- as.character(attrs["condition_type"]),
                   "ShortName" = isSName <<- TRUE,
                   "Organism" = isORG <<- TRUE,
                   "ComplexComponent" = compId <<- as.integer(attrs["idref"]),
                   "Family" = doFamily(attrs),
                   "Interaction" = interID <<- paste(as.character(attrs["id"]),as.character(attrs["interaction_type"]),sep="@"),
                   "Source" = isSource <<- TRUE,
                   "InteractionComponent" =  interCompId <<- paste(as.integer(attrs["idref"]),as.character(attrs["role_type"]),sep="@"),
                   "MoleculeList" = isMol <<- TRUE,
                   "InteractionList" = isInter <<- TRUE,
                   "Pathway" = pathId <<- as.integer(attrs[["id"]]),
                   "PathwayList" = isPath <<- TRUE,
                   "PathwayComponent" = doPathComp(attrs))

        }
        endElement <- function(name,...){
            #print(paste("end", name))
            switch(name,
                   "Molecule" =  writeMol(),
                   "Interaction" =  writeInter(),
                   "InteractionComponent" = writeInterComp(),
                   "Pathway" = writePath(),
                   "ComplexComponent" = writeComp())
        }

        list(startElement = startElement, endElement = endElement,
             text = text)
    }
    xmlEventParse(sourceFile, handlers = getHandler())

    return(list(molecule = molList, interaction = interList,
                pathway = pathList))
}

#getCompnent <- function(compVect){
#    if(is.na(compVect[1])){
#        return(NA)
#    }else{
#        # unique ignores names and has to do this way
#        temp <- unique(paste(compVect, names(compVect), sep = "@"))
#        return(twoStepSplit(paste(temp, collapse = ";"), asNumeric = TRUE))
#    }
#}
