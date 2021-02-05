# Parses the descXXXX file from MeSH, where XXXX are numerics for
# year.
# mesh - the name for the descXXXX file downloaded from MeSH and
#        stored locally.
#

MeSHParser <- function(mesh){
    treeEnv <- new.env(hash = TRUE, parent = NULL)
    scopeEnv <- new.env(hash = TRUE, parent = NULL)
    qualifyEnv <- new.env(hash = TRUE, parent = NULL)
    conceptEnv <- new.env(hash = TRUE, parent = NULL)
    termEnv <- new.env(hash = TRUE, parent = NULL)
    dHeadingEnv <- new.env(hash = TRUE, parent = NULL)
    cHeadingEnv <- new.env(hash = TRUE, parent = NULL)
    tHeadingEnv <- new.env(hash = TRUE, parent = NULL)
    qHeadingEnv <- new.env(hash = TRUE, parent = NULL)

    dHeading <- NULL
    cHeading <- NULL
    tHeading <- NULL
    qHeading <- NULL
    dUid <- NULL
    cUid <- NULL
    tUid <- NULL
    qUid <- NULL
    tree <- NULL
    scope <- NULL
    text <- NULL
    headings <- NULL
    whichOne <- setVars()

    init <- function(){
        dUid <<- NULL
        cUid <<- NULL
        tUid <<- NULL
        qUid <<- NULL
        tree <<- NULL
        scope <<- NULL
    }
    doHeading <- function(){
        multiassign(dHeading[,1], dHeading[,2], dHeadingEnv)
        multiassign(cHeading[,1], cHeading[,2], cHeadingEnv)
        multiassign(tHeading[,1], tHeading[,2], tHeadingEnv)
        multiassign(qHeading[,1], qHeading[,2], qHeadingEnv)
    }
    getHandler <- function(){
        doString <- function(){
            varName <- whichOne[whichOne]
            if(length(varName) > 0){
                print(text)
                switch(names(varName),
                   "dUid" = dHeading <<- rbind(dHeading,
                                            c(dUid[length(dUid)], text)),
                   "cUid" = cHeading <<- unique(rbind(cHeading,
                                            c(cUid[length(cUid)], text))),
                   "tUid" = tHeading <<- unique(rbind(tHeading,
                                            c(tUid[length(tUid)], text))),
                   "qUid" = qHeading <<- unique(rbind(qHeading,
                                            c(qUid[length(qUid)], text))))
                whichOne <<- setVars()
                text <<- NULL
            }
        }
        endElement <- function(name, attrs){
            switch(name,
                   "DescriptorUI" = dUid <<- text,
                   "DescriptorName" = whichOne["dUid"] <<- FALSE,
                   "QualifierUI" = qUid <<- c(qUid, text),
                   "QualifierName" = whichOne["qUid"] <<- FALSE,
                   "AllowableQualifiersList" = assign(dUid, qUid, qualifyEnv),
                   "TreeNumber" = tree <<- c(tree, text),
                   "TreeNumberList" = assign(dUid, tree, treeEnv),
                   "ConceptName" =  whichOne["cUid"] <<- FALSE,
                   "ConceptUI" = cUid <<- c(cUid, text),
                   "Concept" = assign(dUid, cUid, conceptEnv),
                   "ScopeNote" = assign(dUid, text, scopeEnv),
                   "TermUI" = tUid <<- c(tUid, text),
                   "Term" = whichOne["tUid"] <<- FALSE,
                   "TermList" = assign(cUid,tUid, termEnv),
                   "String" = doString(),
                   "DescriptorRecordSet" = doHeading())
        }
        startElement <- function(name, ...){
            switch(name,
                   "DescriptorRecord" = init(),
                   "AllowableQualifierList" = qUid <<- NULL,
                   "Concept" = cUid <<- NULL,
                   "TermList" = tUid <<- NULL,
                   "Term" = whichOne["tUid"] <<- TRUE,
                   "DescriptorName" = whichOne["dUid"] <<- TRUE,
                   "ConceptName" = whichOne["cUid"] <<- TRUE,
                   "QualifierName" = whichOne["qUid"] <<- TRUE)
        }
        text <- function(x){
            if(x != ""){
                text <<- sub("^ ", "", paste(text, x))
            }else{
                text <<- NULL
            }
        }
        list(startElement = startElement, endElement = endElement,
             text = text)
    }
    xmlEventParse(mesh, handlers = getHandler())
    return(list(treenum = treeEnv, scopenote = scopeEnv,
                qualifier = qualifyEnv, concept = conceptEnv,
                term = termEnv, dHeading = dHeadingEnv,
                cHeading = cHeadingEnv, tHeading = tHeadingEnv,
                qHeading = qHeadingEnv))
}

setVars <- function(){
    temp <- rep(FALSE, 4)
    names(temp) <- c("dUid", "cUid", "tUid", "qUid")
    return(temp)
}







