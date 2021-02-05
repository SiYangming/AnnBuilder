#mapLLNGB <- function(organism, pkgName, pkgPath,
#                     ugUrl = getSrcUrl("ug", organism),
#                     egUrl = paste(getSrcUrl("eg"), "gene2accession.gz"),
#                     fromWeb = TRUE){
#
#    repList <- getRepList4Perl(organism, ugUrl, egUrl, fromWeb)
#    llNGB <- getLLNGBMap(repList, "ll2gb")
#    gbNLL <- getLLNGBMap(repList, "gb2ll")
#    saveMat(llNGB, fun = splitEntry, pkgName = pkgName,
#                pkgPath = pkgPath, envName = "LL2ACCNUM")
#    saveMat(gbNLL, fun = function(x) splitEntry(x, asNumeric = TRUE),
#            pkgName = pkgName, pkgPath = pkgPath, envName = "ACCNUM2LL")
#
#    if(fromWeb){
#        unlink(c(repList[["UGFILE"]], repList[["LLFILE"]]))
#    }
#}

mapLLNGB <- function(organism, pkgName, pkgPath,
                     ugUrl = getSrcUrl("ug", organism),
                     egUrl = paste(getSrcUrl("eg"), "gene2accession.gz"),
                     fromWeb = TRUE){

    repList <- getRepList4Perl(organism, ugUrl, egUrl, fromWeb=FALSE)
    
    eg = read.table(repList$LLFILE, na.strings="-", sep="\t", quote="", comment.char="#", header=FALSE)
    
    eg_out = eg[,c(1,2,4,6)]
    eg_out = eg_out[eg_out[,1]==getTaxid(organism),]
    col3 = lapply(as.list(eg_out[,3]), function(x) strsplit(x,"\\.")[[1]][1])
    col4 = lapply(as.list(eg_out[,4]), function(x) strsplit(x,"\\.")[[1]][1])
    dat = rbind(cbind(eg_out[,2],col3), cbind(eg_out[,2],col4))
    dat = dat[!is.na(dat[,2]),]
    
    gbNLL <- split(as.numeric(dat[,1]), as.character(dat[,2]))
    llNGB  <- split(as.character(dat[,2]), as.character(dat[,1]))

    saveList(llNGB, pkgName = pkgName, pkgPath = pkgPath, envName = "LL2ACCNUM")
    saveList(gbNLL, pkgName = pkgName, pkgPath = pkgPath, envName = "ACCNUM2LL")

    if(fromWeb){
        unlink(c(repList[["UGFILE"]], repList[["LLFILE"]]))
    }
}

getTaxid <- function(organism)
{
    switch(organism,
    	"Homo sapiens" = return("9606"),
    	"Rattus norvegicus" = return("10116"),
    	"Mus musculus"= return("10090"),
    	stop("Unknown entry for organism"))

}

mapUGNGB <- function(organism, pkgName, pkgPath,
                     ugUrl = getSrcUrl("ug", organism),
                     llUrl = getSrcUrl("ll"), fromWeb = TRUE){

    repList <- getRepList4Perl(organism, ugUrl, llUrl, fromWeb)
    ugNGB <- getLLNGBMap(repList, "ug2gb")
    gbNUG <- getLLNGBMap(repList, "gb2ug")
    saveMat(ugNGB, fun = splitEntry, pkgName = pkgName,
                pkgPath = pkgPath, envName = "UG2ACCNUM")
    saveMat(gbNUG, fun = splitEntry,
            pkgName = pkgName, pkgPath = pkgPath, envName = "ACCNUM2UG")

    if(fromWeb){
        unlink(c(repList[["UGFILE"]], repList[["LLFILE"]]))
    }
}


getRepList4Perl <- function(organism, ugUrl = getSrcUrl("ug", organism),
                     llUrl = getSrcUrl("ll"), fromWeb = TRUE){

    repList <- list()
    if(fromWeb){
        tempUG <- loadFromUrl(ugUrl)
        tempLL <- loadFromUrl(llUrl)
    }else{
        tempUG <- ugUrl
        tempLL <- llUrl
    }
    if(.Platform$OS.type == "unix"){
        perlBin <- system("which perl", intern = TRUE)
        if(length(perlBin) > 1){
            stop("Perl is not available!")
        }
        repList[["PERLLOC"]] <- paste("#!", perlBin, sep = "")
    }else if(OS == "windows"){
        repList[["PERLLOC"]] <- "#!/usr/bin/perl -w"
    }
    orgCode <- getOrgNameNCode()
    repList[["UGFILE"]] <- tempUG
    repList[["LLFILE"]] <- tempLL
    #repList[["ORGANISM"]] <- getOrgName(organism, "scientific")
    repList[["ORGANISM"]] <- organism
    repList[["ORGCODE"]] <- names(orgCode[orgCode == organism])

    return(repList)
}

# Does two way mappings between LocusLink ids and GenBank Accession
# numbers. Source file used is from UniGene
# (ftp://ftp.ncbi.nih.gov/repository/UniGene).
# ug - a UG object with correct url, parser, ... set.

getLLNGBMap <- function(repList, what = "ll2gb"){
    tempPerl <- paste(tempfile("ugNllParser"), "pl", sep=".")
    repList[["OUTFILE"]] <- tempfile("llNgb")
    switch(toupper(what),
           LL2GB = parser <- "ll2GBPerl.temp",
           GB2LL = parser <- "gb2LLPerl.temp",
           UG2GB = parser <- "ug2GBPerl.temp",
           GB2UG = parser <- "gb2UGPerl.temp",
           stop("Unknown mapping type"))

    copySubstitute(file.path(.path.package("AnnBuilder"),
                            "templates", parser), tempPerl, repList, "#")
    .callPerl(tempPerl, .Platform$OS.type)

    parsed <- matrix(scan(repList[["OUTFILE"]], what = "character",
                           sep = "\t", quote = "", quiet = TRUE,
                           strip.white = TRUE, comment.char = ""),
                           ncol = 2, byrow = TRUE)
    unlink(c(tempPerl, repList[["OUTFILE"]]))
    return(parsed)
}
