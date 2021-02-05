# Get a matrix with Affymetrix
yeastAnn <- function(base = "", yGenoUrl,
    yGenoNames = c("literature_curation/gene_literature.tab",
      "chromosomal_feature/SGD_features.tab",
      "literature_curation/gene_association.sgd.gz"), toKeep =
                     list(c(6, 1), c(1, 5, 9, 10, 12, 16, 6), c(2, 5, 7)),
                     colNames = list(c("sgdid", "pmid"),
                     c("sgdid", "genename", "chr", "chrloc", "chrori",
                       "description", "alias"),
                     c("sgdid", "go")),
                     seps = c("\t", "\t", "\t"),
                     by = "sgdid"){
    if (missing(yGenoUrl))
      yGenoUrl <- getSrcUrl("YG")

    # Each file should have matching parsing parameters
    if(any(c(length(toKeep), length(seps), length(colNames)) !=
           length(yGenoNames))){
        stop("Lengths of fileNames, toKeep, and seps have to be the same")
    }

    temp1 <- getProbe2SGD(base, yGenoUrl = yGenoUrl)
    # Creat a yeastGeno object
    ygeno <- YG(srcUrl = yGenoUrl)
    # merge temp1 with each of the annotation data defined by fileNames
    for(i in 1:length(yGenoNames)){
        temp <- readData(ygeno, yGenoNames[i], toKeep[[i]], seps[i])
        if(length(grep(".*gene_association.sgd", yGenoNames[i])) > 0){
            temp <- mergeRowByKey(cbind(temp[,1],
                                 paste(temp[,2], temp[,3], sep = "@")))
            colnames(temp) <- colNames[[i]]
             # Remove the leading and ending space
            dataCName <- colnames(temp)
            temp <- t(apply(temp, 1, function(x) gsub("^ | $", "", x)))
            colnames(temp) <- dataCName
            #temp[,"go"] <- formatGO(temp[,"go"], temp[,"evi"])
            #temp <- temp[, setdiff(colNames[[i]], "evi")]

        }else if(length(grep(".*SGD.*tab$", yGenoNames[i])) > 0){
            colnames(temp) <- colNames[[i]]
            # Remove the leading and ending space
            dataCName <- colnames(temp)
            temp <- t(apply(temp, 1, function(x) gsub("^ | $", "", x)))
            colnames(temp) <- dataCName
            temp[,"chrloc"] <- formatChrLoc(temp[,"chr"], temp[,"chrloc"],
                                            temp[,"chrori"])
            temp[, "alias"] <- gsub("\\|", ";", temp[, "alias"])
            ##temp[, "alias"] <- sapply(temp[, "alias"], function(x) ifelse(x=="", NA, x), USE.NAMES=FALSE)
            ##temp[, "genename"] <- sapply(temp[, "genename"], function(x) ifelse(x=="", NA, x), USE.NAMES=FALSE)
            temp <- temp[, setdiff(colNames[[i]], "chrori")]
        }else{
            colnames(temp) <- colNames[[i]]
             # Remove the leading and ending space
            dataCName <- colnames(temp)
            temp <- t(apply(temp, 1, function(x) gsub("^ | $", "", x)))
            colnames(temp) <- dataCName
        }
        # Only merge ones mapped to SGD ids
        temp1 <- merge(temp1, mergeRowByKey(temp), by = by, all.x = TRUE)
#        temp1 <- (temp1[!duplicated(temp1[,by]) |
#                                         !duplicated(temp1[,"probe"]),])
    }

    ##improve the genename / alias
    temp <- readData(ygeno, "gene_registry/registry.genenames.tab", c(1,2,7), "\t")
    colnames(temp) <- c("genename2", "alias2", "sgdid")
    temp[, "alias2"] <- gsub("\\|", ";", as.vector(temp[, "alias2"]))
    temp1 <- merge(temp1, mergeRowByKey(temp), by = by, all.x = TRUE)
    ## an ugly patch for genename (Oct, 2006)
    temp1[, "genename"] <- as.character(temp1[, "genename"])
    temp1[, "genename2"] <- as.character(temp1[, "genename2"])
    naGeneNameIndex <- which(temp1[,"genename"]=="")
    temp1[naGeneNameIndex, "genename"] <- temp1[naGeneNameIndex,"genename2"]

    temp1[, "alias"] <- apply(temp1[, c("alias", "alias2")], 1,
                              function(x){
                                  xx <- union(unlist(strsplit(x[1], ";")), unlist(strsplit(x[2], ";")))
                                  if(all(is.na(xx))){
                                      return(NA)
                                  }else{
                                      xx <- xx[!is.na(xx)]
                                      return(paste(xx, sep="", collapse=";"))
                                  }
                              })
    ## Return data with ORF removed
    return(temp1[, -match(c(by, "genename2", "alias2"), colnames(temp1))])
}
# Gets mappings between ORF to sgdid
getProbe2SGD <- function(probe2ORF = "", yGenoUrl,
                     fileName = "literature_curation/orf_geneontology.tab",
                     toKeep = c(1, 7), colNames = c("orf", "sgdid"),
                     sep = "\t", by = "orf"){
    yGeno <- YG(srcUrl = yGenoUrl)
    temp <- readData(yGeno, fileName, toKeep, sep)
    colnames(temp) <- colNames
    if(probe2ORF == ""){
        colnames(temp) <- c("probe", "sgdid")
        return(temp)
    }else{
        probes <- read.table(probe2ORF, sep = "\t", header = FALSE,
		stringsAsFactors=FALSE,
            strip.white = TRUE)
        colnames(probes) <- c("probe", by)
        merged <- merge(probes, temp, by = by, all.x = TRUE)
        # Remove the duplicates. Keep the first one
        merged <- merged[!duplicated(merged[,1]),]
        return(merged)
    }
}

# A function to extract yeast annotation data from
# "ftp://genome-ftp.stanford.edu/pub/yeast/data_download".
#
# Data will be extracted from four files each with an column for
# source specific ids linking the data together.
#
# Copyright 2003, J. Zhang, all rights reserved.
#
procYeastGeno <- function(baseURL, fileName, toKeep, colNames, seps = "\t"){

    # Parse each file and merge based on source specified ids
    merged <- NULL
    # Creat a yeastGeno object
    ygeno <- YG(srcUrl = baseURL)
    for(i in 1:length(fileName)){
        temp <- readData(ygeno, fileName[i], toKeep[[i]], seps[i])
        # The first row is for column name. Remove it
        temp <- temp[-1,]
        colnames(temp) <- colNames[[i]]
        if(is.null(merged)){
            merged <- temp
        }else{
            # Only merge the ones mapped to GenBank accessin number
            merged <- merge(merged, temp, by = "id", all = TRUE)
        }
    }
    # Drop the source specific id
    return(merged)
}
# GO ids from yeast data source do not have the leading GO and 0s. Put
# them back
formatGO <- function(gos, evis){
    add0 <- function(go){
        num <- 7 - nchar(go)
        return(paste("GO:", paste(rep("0", num), sep = "",
                                  collapse = ""), go, sep = "", collapse = ""))
    }
    gos <- sapply(gos, add0)
    return(paste(gos, evis, sep = "@"))
}
# Name chromosomal locations by chromosome number and add + for W and
# - for C strand
formatChrLoc <- function(chr, chrloc, chrori){
    chrori <- sapply(chrori, function(x) ifelse(x == "W", "+", "-"))
    chrloc <- paste(chrori, chrloc, "@", chr, sep = "")
    return(chrloc)
}

# Get the data that map Affymetrix probe ids to SGD_ORF identifiers
# that in trun can be mapped to annotation data using the source from
# yeast genomic web site.
getGEOYeast <- function(GEOAccNum, GEOUrl, geoCols = c(1, 8), yGenoUrl) {
    geo <- GEO(GEOUrl)
    geoData <- readData(geo, GEOAccNum)
    # Read the file with mappigs between Affymetrix probe ids and
    # SGD_ORF identifiers
    geoData <- geoData[,geoCols]
    colnames(geoData) <- c("probe", "orf")
    return(geoData)
}

getYGExons <- function(srcUrl, yGenoName = "chromosomal_feature/intron_exon.tab",
                       sep = "\t") {
    ygeno <- YG(srcUrl = srcUrl)
    temp <- getProbe2SGD("", yGenoUrl = srcUrl)
    exons <- readData(ygeno, yGenoName, c(3, 5, 9, 10), sep)
    colnames(exons) <- c("sgdid", "chrom", "start", "end")
    exons <-as.matrix( merge(temp, exons, by = "sgdid", all.x = TRUE))
    # Remove NA
    exons <- exons[!is.na(exons[, "chrom"]), ]
    # Remove white spaces
    for(i in c("chrom", "start", "end")){
        exons[, i] <- gsub(" ", "", exons[, i])
    }

    return(exons[, c("probe", "chrom", "start", "end")])
}







