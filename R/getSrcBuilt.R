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
           "PFAM" = return(getPFAMBuilt()),
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

getUGBuilt <- function(organism){
  url <- .srcUrls("UG")

    #switch(toupper(organism),
    #       "HUMAN" = infoUrl <- paste(url, "/Hs.info", sep = ""),
    #       "MOUSE" = infoUrl <- paste(url, "/Mm.info", sep = ""),
    #       "RAT" = infoUrl <- paste(url, "/Rn.info", sep = ""),
    #       stop(paste("Organism", organism, "is not supported.")))

    infoUrl <- paste(url, "/", gsub(" ", "_", organism), "/",
                     getUGShortName(organism), ".info", sep = "")

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

    ## This part is removed by Ting-Yuan because the "built" is included in the url
    ##switch(toupper(organism),
    ##       "MUS MUSCULUS" =  built <- gsub(".*goldenPath/mm(.*)/database.*",
    ##                       "\\1", url),
    ##       "HOMO SAPIENS" = built <-
    ##                      gsub(".*goldenPath/(.*)/database.*", "\\1", url),
    ##       "RATTUS NORVEGICUS" = built <-
    ##                      gsub(".*goldenPath/rn(.*)/database.*", "\\1", url),
    ##       "DANIO RERIO" = built <-
    ##                 gsub("goldenPath/danRer(.*)/database.*", "\\1" , url),
    ##       return("N/A"))
    ##       #stop(paste("Organism", organism, "is not supported")))
    ##
    ##switch(toupper(organism),
    ##       "HOMO SAPIENS" = built <- gsub(".*goldenPath/(.*)/database.*",
    ##                       "\\1", url),
    ##       "MUS MUSCULUS" = key <- "goldenPath/mm.*Annotation database",
    ##       "RATTUS NORVEGICUS" = key <- "goldenPath/rn.*Annotation database",
    ##       "DANIO RERIO" = key <- "goldenPath/danRer.*Annotation database",
    ##       key <- NA)

    built <- "No build info available."
    
    if(is.na(built) || is.null(built) || built == ""){
        warning("Built for UCSC is not valid!")
        return("N/A")
    }else{
        return(built)
    }
}

getGOBuilt <- function(){
  dateStrings <- try(readLines(dirname(getSrcUrl("GO"))), TRUE)
  buildInfo <- try(unlist(strsplit(dateStrings[grep("-termdb.rdf-xml.gz",
                                           dateStrings)], "</a>")), TRUE)
  build <- try(gsub(" *([0-9a-zA-Z-]+) .*", "\\1",
                    buildInfo[length(buildInfo)]), TRUE)
  if(any(c(class(dateStrings), class(buildInfo), class(build))
         == "try-error")){
    return("Build Information not available")
  }else{
    return(build)
  }
    #goUrl <- getSrcUrl("GO")
    #datePart <- gsub("^.*go_([0-9]+)-.*$", "\\1", goUrl)
    #datePart
}

getKEGGBuilt <- function(
               url = "http://www.genome.jp/kegg/docs/relnote.html"){
    options(show.error.messages = FALSE)
    got <- try(readURL(url))
    options(show.error.messages = TRUE)
    if(inherits(got, "try-error")){
        warning(paste("Can't read", url))
        return("N/A")
    }
    aLine <- got[grep("^Release", got)][1]
    built <- gsub(".*(Release [0-9.]* \\(.*\\)) .*", "\\1", aLine)
    if(is.na(built) || is.null(built) || built == "" || length(built)==0 ){
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

getPFAMBuilt <- function(){
    return("PFAM built data not available")
}
