# Functions that return the built date for a given public source.
#
# Copyright 2002, J. Zhang, all rights reserved
#

getSrcBuilt <- function(src = "LL", organism = "Homo sapiens"){
    switch(toupper(src),
           "LL" = return(getLLBuilt()),
           "UG" = return(getUGBuilt(organism)),
           "GP" = return(getUCSCBuilt(organism)),
           "GO" = return(getGOBuilt()),
           "KEGG" = return(getKEGGBuilt()),
           "YG" = return(getYGBuilt()),
           "HG" = return(getHGBuilt()),
           "REFSEQ" = return(getRefSeqBuilt(organism)),
           "EG" = return(getEGBuilt()),
           return("Build information unavailable"))
}

getEGBuilt <- function(){
  return(paste("Source data downloaded from Entrez Gene on", date()))
}

getLLBuilt <- function(
        url = "http://www.ncbi.nlm.nih.gov/LocusLink/statistics.html"){


    ## <FIXME>
    ## got <- readURL(url)
    ## built <- gsub(".*LocusLink Statistics dated (.*)",
    ##               "\\1", got[grep("LocusLink Statistics dated", got)])
    ## XXX: The webpage with data data is not currently available. Instead
    ##      we will use the run date.
    built <- date()
    ## </FIXME>
    if(is.na(built) || is.null(built) || built == ""){
        warning("Built for LL is not valid!")
        return("N/A")
    }else{
        return(built)
    }
}

getUGBuilt <- function(organism,
                       url = "ftp://ftp.ncbi.nih.gov/repository/UniGene"){

    #switch(toupper(organism),
    #       "HUMAN" = infoUrl <- paste(url, "/Hs.info", sep = ""),
    #       "MOUSE" = infoUrl <- paste(url, "/Mm.info", sep = ""),
    #       "RAT" = infoUrl <- paste(url, "/Rn.info", sep = ""),
    #       stop(paste("Organism", organism, "is not supported.")))

    infoUrl <- paste(url, "/", getUGShortName(organism),
                 ".info", sep = "")

    got <- readURL(infoUrl)
    built <- gsub(".*(Build #[0-9]*).*", "\\1",
                  got[grep("Build #", got)[1]])
    if(is.na(built) || is.null(built) || built == ""){
        warning("Built for UG is not valid!")
        return("N/A")
    }else{
        return(built)
    }
}

getUCSCBuilt <- function(organism){
    url <- getSrcUrl(src = "GP", organism = organism)
    switch(toupper(organism),
           "MUS MUSCULUS" =  built <- gsub(".*goldenPath/mm(.*)/database.*",
                           "\\1", url),
           "HOMO SAPIENS" = built <-
                          gsub(".*goldenPath/(.*)/database.*", "\\1", url),
           "RATTUS NORVEGICUS" = built <-
                          gsub(".*goldenPath/rn(.*)/database.*", "\\1", url),
           "DANIO RERIO" = built <-
                     gsub("goldenPath/danRer(.*)/database.*", "\\1" , url),
           return("N/A"))
           #stop(paste("Organism", organism, "is not supported")))

    switch(toupper(organism),
           "HOMO SAPIENS" = built <- gsub(".*goldenPath/(.*)/database.*",
                           "\\1", url),
           "MUS MUSCULUS" = key <- "goldenPath/mm.*Annotation database",
           "RATTUS NORVEGICUS" = key <- "goldenPath/rn.*Annotation database",
           "DANIO RERIO" = key <- "goldenPath/danRer.*Annotation database",
           key <- NA)

    if(is.na(built) || is.null(built) || built == ""){
        warning("Built for UCSC is not valid!")
        return("N/A")
    }else{
        return(built)
    }
}

getGOBuilt <- function(){
    goUrl <- getSrcUrl("GO")
    datePart <- gsub("^.*go_([0-9]+)-.*$", "\\1", goUrl)
    datePart
}

getKEGGBuilt <- function(url = "http://www.genome.ad.jp/kegg/kegg2.html"){
    options(show.error.messages = FALSE)
    got <- try(readURL(url))
    options(show.error.messages = TRUE)
    if(inherits(got, "try-error")){
        warning(paste("Can't read", url))
        return("N/A")
    }
    aLine <- got[grep("<i>KEGG Release", got)]
    built <- gsub(".*(Release [0-9.]* \\(.*\\)) .*", "\\1", aLine)
    if(is.na(built) || is.null(built) || built == ""){
        warning("Built for KEGG is not valid!")
        return("N/A")
    }else{
        return(built)
    }
}

getYGBuilt <- function(){
    return(paste("Yeast Genome data are built at various time ",
                 "intervals. Sources used were downloaded ", date(), "\n"))
}

readURL <- function(url){
    con <- url(url)
    options(show.error.messages = FALSE)
    temp <- try(readLines(con))
    close(con)
    options(show.error.messages = TRUE)
    if(!inherits(temp, "try-error")){
        return(temp)
    }else{
        stop(paste("Can't read from url:", url))
    }
}

getRefSeqBuilt <- function(organism){
    switch(tolower(organism),
           "human" = return(getRefBuilt4HS()),
           "yeast" = return("Build for yeast not available"),
           stop(paste("Organism", organism, "not supported")))
}

getHGBuilt <- function(){
    return("HomoloGene built date not available")
}

getRefBuilt4HS <- function(){
    url <- getSrcUrl("REFSEQ", "human")
    url <- paste(url, "README_CURRENT_BUILD", sep = "")
    text <- readLines(url)
    return(text[grep("NCBI Build Number:.*", text)])
}











