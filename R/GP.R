# Constructs an object of GP for Golden Path data
#

GP <- function(srcUrl = getSrcUrl("GP"), organism = "human",
               parser = "", baseFile = "", built = getSrcBuilt("gp",
               organism = organism), fromWeb = TRUE){

    new("GP", srcUrl = srcUrl, organism = organism, parser = parser,
        baseFile = baseFile, built = built, fromWeb = TRUE)
}

getChroLocation <- function(srcUrl, exten = gpLinkNGene(), sep = "\t",
                            fromWeb = TRUE, raw = FALSE){
    if (is.na(srcUrl))
	stop("No data from USCS Genome Browser (Golden Path).")	
    if(fromWeb){
        # Tow source data files are needed.
        linkData <- getGPData(paste(srcUrl, exten["link"], sep = ""))
        linkData <- linkData[, c(3, 7)]
        locatData <- getGPData(paste(srcUrl, exten["gene"], sep = ""))
    }else{
        linkData <- read.table(file.path(srcUrl, srcUrl["link"]),
				stringsAsFactors=FALSE,
                               sep = sep, header = FALSE, as.is = TRUE)
        linkData <- linkData[, c(3, 7)]
        locatData <- read.table(file.path(srcUrl, srcUrl["gene"]),
                                sep = sep, header = FALSE,
                                strip.white = TRUE,
				stringsAsFactors=FALSE,
                                as.is = TRUE)
    }
    locatCol <- ncol(locatData)
    if(locatCol==10)
		locatData <- locatData[,c(1:5)]
    else if (locatCol==16)
		locatData <- locatData[,c(2:6)]
    else
		stop("The file format of ", locatData, " cannot be recognized by getChroLocation(). Expect to see 10 or 16 columns, but see ", locatCol, " instead.")
   
    colnames(linkData) <- c("ID", "ENTREZID")
    # Remove "chr" proceeding chromosome number
    locatData[, 2] <- gsub("^chr", "\\", locatData[, 2])
    if(raw){
        colnames(locatData) <- c("ID", "CHR", "STRAND", "START", "END")
        merged <- as.matrix(merge(linkData, locatData, by = "ID",
                                  all.x = TRUE))
        merged <- merged[, 2:6]
        # Remove <NA>s that are introduced by droping the ID column
        merged <- merged[!is.na(merged[,1]), ]
        return(merged)
    }
    # for - strand, end locations should be strat locations
    #locatData <- rbind(locatData[locatData[,3] == "+", c(1:4)],
    #                   locatData[locatData[,3] == "-", c(1, 2, 3, 5)])
    locatData <- locatData[, c(1:4)]
    # Put sign to locations
    locatData[, 4] <- paste(locatData[,3], locatData[,4], sep = "")
    locatData <- cbind(locatData[, 1], paste(locatData[, 4], "@",
                 locatData[,2], sep = ""))
    colnames(locatData) <- c("ID", "CHRLOC")
    merged <- as.matrix(merge(linkData, locatData, by = "ID",
                              all.x = TRUE)[, 2:3])
    # Remove <NA>s that are introduced by droping the ID column
    merged <- merged[!is.na(merged[,1]), ]

    ## clean up
    if (fromWeb)
      unlink(c(linkData, locatData))
    return(mergeRowByKey(merged))
}

getGPData <- function(srcUrl, sep = "\t"){
    temp <- loadFromUrl(srcUrl)
    oneRecord <- scan(temp, what="character", nmax=1, sep="\n")
    ncol <- length(unlist(strsplit(oneRecord, "\t")))
    gpData <- matrix(scan(temp, what = "", sep = sep, quote = "",
                          strip.white = TRUE, quiet = TRUE), ncol = ncol,
                     byrow = TRUE)
    unlink(temp)
    return(gpData)
}

gpLinkNGene <- function(test = FALSE, fromWeb = TRUE){
    if(test){
        return(c(link = "Tlink.txt", gene = "TGene.txt"))
    }else{
        if(fromWeb){
            return(c(link = "refLink.txt.gz", gene = "refGene.txt.gz"))
        }else{
            return(c(link = "refLink.txt", gene = "refGene.txt"))
        }
    }
}

# No in use at this time
#gpParser <- function(){
#    c(link = file.path(.path.package("pubRepo"), "data", "gpLinkParser"),
#      gene = file.path(.path.package("pubRepo"), "data", "gpGeneParser"))
#}

getCytoLoc <- function(organism, srcUrl = paste(getSrcUrl("gp",
          organism), "/",
          "http://www.genome.ucsc.edu/goldenPath/mm4/database/" )){
    cytoBand <- list()
    cytoLoc <- getGPData(srcUrl)[,1:4]
    cytoLoc[,1] <- gsub("chr", "", cytoLoc[,1])
    for(i in unique(cytoLoc[,1])){
        cytoBand[[i]] <- getCytoList(cytoLoc[cytoLoc[,1] == i,
                                             2:ncol(cytoLoc)])
    }
    return(cytoBand)
}

getCytoList <- function(data){
    getStartNEnd <- function(band){
        return(c("Start" = as.numeric(data[data[,3] == band,1]),
               "End" = as.numeric(data[data[,3] == band,2])))
    }
    if(is.null(data)){
        return(NULL)
    }else{
        if(is.null(nrow(data))){
            data <- matrix(data, ncol = 3)
        }
    }

    tempList <- lapply(data[,3], getStartNEnd)
    names(tempList) <- data[,3]
    return(tempList)
}

