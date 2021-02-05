GO <- function(srcUrl = getSrcUrl("GO", xml = TRUE), parser = "",
               baseFile = "", built, fromWeb = TRUE){

    new("GO", srcUrl = srcUrl, parser = parser, baseFile = baseFile,
        built = ifelse(missing(built), getSrcBuilt("GO"), built),
        fromWeb = fromWeb)
}


GOXMLParser <- function(fileName) 
{
    goObsoleteEnv = new.env(hash=TRUE, parent=NULL)
    goTermEnv = new.env(hash=TRUE, parent=NULL)
    cbpEnv = new.env(hash=TRUE, parent=NULL)
    cmfEnv = new.env(hash=TRUE, parent=NULL)
    cccEnv = new.env(hash=TRUE, parent=NULL)
    bpEnv = new.env(hash=TRUE, parent=NULL)
    ccEnv = new.env(hash=TRUE, parent=NULL)
    mfEnv = new.env(hash=TRUE, parent=NULL)
    
    goXmlHandlers <- function() {
        isGOID <- FALSE
        isName <- FALSE
        isCat <- FALSE
        isAltID <- FALSE
        isType <- FALSE
        isTO <- FALSE
        isa <- FALSE
        isObsolete <- FALSE
        isSyn <- FALSE
        isDef <- FALSE
        
        GOID <- NA
        parents <- character(0)
        goname <- NA
        gocat <- NA
        alt_id <- character(0)
        datatype <- character(0)
        related_id <- character(0)
        obsolete <- NA
        definition <- NA
        synonym <- NA

        setDefault <- function()
          {
              GOID <<- NA
              parents <<- character(0)
              goname <<- NA
              gocat <<- NA
              definition <<- NA
              synonym <<- NA
              alt_id <<- character(0)
              datatype <<- character(0)
              related_id <<- character(0) 
              isObsolete <<- FALSE
          }
        
        appendParent <- function(attrs, assoType){
            parent_id = attrs
            names(parent_id) = assoType
            parents <<- c(parents, parent_id)
        }
        
        
        startElement <- function(name, attrs){
            switch(name,
                   "id" = isGOID <<- TRUE,
                   "name" = isName <<- TRUE,
                   "namespace" = isCat <<- TRUE,
                   "defstr" = isDef <<- TRUE,
                   "synonym_text" = isSyn <<- TRUE,
                   "alt_id" = isAltID <<- TRUE,
                   "type" = isType <<- TRUE,
                   "to" = isTO <<- TRUE,
                   "is_obsolete" = isObsolete <<- TRUE,
                   "is_a" = isa <<- TRUE)
        }

        endElement <- function(name, ...){
            if (name == "term") {
                parents = getRelation(datatype,related_id)
                if (!isObsolete) {
                    ## Create a GOTerms object from annotate
                    if (length(alt_id) > 0) {
                        for (id in alt_id) {
                            goTerm = GOTerms(
                              GOId=as.character(id),
                              term=as.character(goname[2]),
                              ontology=as.character(gocat), 
                              synonym=as.character(synonym),
                              secondary=character(0),
                              definition=as.character(definition))
                            assign(id, goTerm, goTermEnv)
                        }
                    }
                    goTerm = GOTerms(
                      GOId=as.character(GOID),
                      term=as.character(goname[2]),
                      ontology=as.character(gocat), 
                      synonym=as.character(synonym),
                      secondary=as.character(alt_id),
                      definition=as.character(definition))
                    assign(GOID, goTerm, goTermEnv)
                    if (gocat == "BP") {
                        if (length(alt_id) > 0) {
                            for (id in alt_id) {
                                assign(id, parents, bpEnv)
                            }
                        }
                        assign(GOID, parents, bpEnv)
                    }
                    else if (gocat == "MF") {
                        if(length(alt_id) > 0) {
                            for(id in alt_id) {
                                assign(id, parents, mfEnv)
                            }
                        }
                        assign(GOID, parents, mfEnv)
                    }
                    else if (gocat == "CC") {
                        if (length(alt_id) > 0) {
                            for (id in alt_id) {
                                assign(id, parents, ccEnv)
                            }
                        }
                        assign(GOID, parents, ccEnv)
                    }
                    else {
                        warning(paste("GOID", sQuote(GOID),"contains an unknown category:",sQuote(gocat)))
                        
                    }
                }
                else {
                    if (length(alt_id) > 0) {
                        for (id in alt_id) {
                            goTerm = GOTerms(
                              GOId=as.character(id),
                              term=as.character(goname[2]),
                              ontology=as.character(gocat),
                              definition=as.character(definition))
                            assign(id, goTerm, goObsoleteEnv)
                        }
                    }
                    goterm = GOTerms(
                      GOId=as.character(GOID),
                      term=as.character(goname[2]),
                      ontology=as.character(gocat),
                      definition=as.character(definition))
                    assign(GOID, goterm, goObsoleteEnv)
                }
                
                setDefault()
            }
        }
        
        getRelation <- function(datatype, related_id) {
            if (length(datatype) > 0) {
                names(related_id) = datatype 
                parents <<- c(parents, related_id)
                return(parents)
            }
            else {
                return (parents)
            }         
        }
        
        text <- function(theText) {
            if (isGOID) {
                GOID <<- theText
                isGOID <<- FALSE
            } else if (isName) {
                goname <<- c(goname, theText)
                isName <<- FALSE
            }
            else if (isCat) {
                ## normalize ontology name
                ## if unrecognized, return the text encountered
                cat <- switch(theText,
                              biological_process="BP",
                              molecular_function="MF",
                              cellular_component="CC",
                              theText)
                gocat <<- cat
                isCat <<- FALSE
            }
            else if (isAltID) {
                alt_id <<- c(alt_id, theText)
                isAltID <<- FALSE
            }
            else if (isType) {
                datatype <<- c(datatype, theText)
                isType <<- FALSE
            }
            else if (isTO) {
                related_id <<- c(related_id, theText)
                isTO <<- FALSE
            }
            else if (isa) {
                appendParent(theText, "isa")
                isa <<- FALSE
            }
            else if (isDef) {
                definition <<- theText
                isDef <<- FALSE
            }
            else if(isSyn) {
                synonym <<- theText
                isSyn <<- FALSE
            }
        }
        
        return(list(startElement=startElement, endElement=endElement, text=text,
                    getbp=function() bpEnv,
                    getmf=function() mfEnv,
                    getcc=function() ccEnv,
                    getObsolete=function() goObsoleteEnv,
                    getGoterm = function() goTermEnv))
    }
    
    results <- xmlEventParse(fileName, handlers=goXmlHandlers())
    ccEnv = addRoot(ccEnv, "cc")
    bpEnv = addRoot(bpEnv, "bp")
    mfEnv = addRoot(mfEnv, "mf")
    
    cbpEnv <- createGOChildren(results$getbp())
    cmfEnv <- createGOChildren(results$getmf())
    cccEnv <- createGOChildren(results$getcc())
    
    obpEnv = createGoOffspring(cbpEnv, "GO:0008150")
    omfEnv = createGoOffspring(cmfEnv, "GO:0003674")
    occEnv = createGoOffspring(cccEnv, "GO:0005575")
    
    abpEnv = createGoAncestor(obpEnv)
    amfEnv = createGoAncestor(omfEnv)
    accEnv = createGoAncestor(occEnv)

    return(list(TERM=goTermEnv,
                BPPARENTS=bpEnv,
                MFPARENTS=mfEnv,
                CCPARENTS=ccEnv,
                BPCHILDREN=cbpEnv,
                MFCHILDREN=cmfEnv,
                CCCHILDREN=cccEnv,
                OBSOLETE=results$getObsolete(),
                BPOFFSPRING=obpEnv,
                MFOFFSPRING=omfEnv,
                CCOFFSPRING=occEnv,
                BPANCESTOR=abpEnv,
                MFANCESTOR=amfEnv,
                CCANCESTOR=accEnv))
}


addRoot <- function(parentEnv, root) {
    rootId <- switch(root,
                     "bp" = "GO:0008150",
                     "cc" = "GO:0005575",
                     "mf" = "GO:0003674"
                     )
    assign(rootId, "all", parentEnv)
    return(parentEnv)
}


createGOChildren <- function(parentEnv) {
    
    numParents <- eapply(parentEnv, length)
    ngoid <- rep(names(as.list(parentEnv)), numParents)
    children <- split(ngoid, unlist(as.list(parentEnv)))
    childEnv = new.env(parent=NULL)
    l2e(children, childEnv)
    
}


createGoOffspring <- function(childrenEnv, rootID) {
    goVec <- unique(unlist(as.list(childrenEnv)))
    goVecKey = unique(unlist(names(as.list(childrenEnv))))
    isChildren = sapply(goVec, function(x) x %in% goVecKey)
    noChildren = new.env(parent=NULL)
    nochildren_ids = goVec[isChildren==FALSE]
    names(nochildren_ids) = nochildren_ids
    l2e(as.list(nochildren_ids),noChildren) 
    GoOffspring = new.env(parent=NULL)
    
    for (y in 1:length(goVec))
      {
          GoOff = findOffspring(goVec[y], noChildren, childrenEnv)
          if (!goVec[y] %in% GoOff) {
              GoOffspring[[goVec[y]]] = findOffspring(goVec[y], noChildren, childrenEnv)
          }
          else {
              GoOffspring[[goVec[y]]] <- NA
          }
      }
    GoOffspring[["all"]] <- c(unique(unlist(as.list(GoOffspring))), rootID)
    return(GoOffspring)
}


findOffspring <- function(goID, mEnv, childrenEnv) {
    found <- mEnv[[goID]]
    if (is.null(found)) {
        children <- childrenEnv[[goID]]
        thidIdsOffspring <- children
        for (childID in children) {
            offspring <- findOffspring(childID, mEnv, childrenEnv)
            mEnv[[childID]] <- offspring
            thidIdsOffspring <- c(thidIdsOffspring, offspring)
        }
    }
    else {
        thidIdsOffspring <- found
    }
    return(unique(unlist(thidIdsOffspring)))
}


createGoAncestor <- function(offspringEnv) {
    numOffspring <- eapply(offspringEnv, length)
    ngoid <- rep(names(as.list(offspringEnv)), numOffspring)
    ancestor <- split(ngoid, unlist(as.list(offspringEnv)))
    ancestorEnv = new.env(parent=NULL)
    l2e(ancestor, ancestorEnv)
}

