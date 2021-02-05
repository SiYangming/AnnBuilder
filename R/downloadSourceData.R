## downloadSourceData.R
##
## Create a local repository of public domain data files (downloaded from the
## internet) for use in building annotation data packages.
##
## Goal is to download everything needed as input data for the annotation
## package builds so that building many packages does not require redownloading
## of the same files over the internet (get them from the LAN instead).
##


writeSourceUrlConfig <- function(file) {
    pubUrls <- getOption("AnnBuilderPublicDataUrls")
    mirrorUrls <- gsub("^[hft]+tp://(.*)$", "\\1", pubUrls)
    ## XXX: KEGG is a special case
    idx <- which(names(pubUrls) == "KEGG")
    mirrorUrls[idx] <- "kegg/pathways"
    df <- data.frame(name=names(pubUrls), url=mirrorUrls)
    write.table(file=file, df, row.names=FALSE, quote=FALSE)
}


readSourceUrlConfig <- function(file, urlPrefix) {
    df <- read.table(file, colClasses=rep("character", 2), header=TRUE)
    urls = df$url
    urls <- as.list(urls)
    if (!missing(urlPrefix))
      urls <- paste(urlPrefix, urls, sep="/")
    urls <- as.list(urls)
    names(urls) = df$name
    urls
}

    
downloadSourceData <- function() {
    pubDataURLs <- getOption("AnnBuilderPublicDataUrls")
    getPubDataHomoloGene(pubDataURLs[["HG"]])
    getPubDataLocusLink(pubDataURLs[["LL"]])
    getPubDataUniGene(pubDataURLs[["UG"]])
    getPubDataEntrezGene(pubDataURLs[["EG"]])
    getPubDataGoldenPath(pubDataURLs[["GP"]])
    getPubDataGo(pubDataURLs[["GO"]])
    getPubDataYeastGenome(pubDataURLs[["YG"]])
    getPubDataKegg(pubDataURLs[["KEGG"]])
}


wget <- function(url, levels, accepts) {
    if (!missing(accepts)) {
        accepts <- paste(accepts, collapse=",")
        accepts <- paste("--accept", accepts)
    } else {
        accepts <- ""
    }
    if (missing(levels))
      levels <- 5
    levels <- paste("--level", levels)
    globalOpts <- paste("--recursive", "--timestamping", "--retr-symlinks",
                        levels, accepts)
        wgetCmd <- paste("wget", globalOpts, url)
    cat("Running:", wgetCmd, "\n")
    system(wgetCmd)
}


getPubDataHomoloGene <- function(baseUrl) {
    wget(baseUrl, levels=1, accepts="hmlg.ftp")
}


getPubDataLocusLink <- function(baseUrl) {
    wget(baseUrl, levels=1, accepts="LL_tmpl.gz")
}


getPubDataUniGene <- function(baseUrl) {
    wget(baseUrl, levels=1, accepts=".data.gz")
}


getPubDataEntrezGene <- function(baseUrl) {
    wget(baseUrl, levels=1)
}


getPubDataGoldenPath <- function(baseUrl) {
    ## Download Golden Path data to current directory
    ##
    ## Get most current genome data for all available species.  Only the data files
    ## listed in wantedFiles are retrieved.
    ##

    ## This is a list of organisms pulled from the golden path website.  Only
    ## some of these have the refGene.txt.gz and refLink.txt.gz files we are
    ## interested in.
    organisms <- c("Anopheles_gambiae", "Apis_mellifera", "Bos_taurus",
                   "Caenorhabditis_briggsae", "Caenorhabditis_elegans",
                   "Canis_familiaris", "Danio_Rerio", "Drosophila_ananassae",
                   "Drosophila_melanogaster", "Drosophila_mojavensis",
                   "Drosophila_pseudoobscura", "Drosophila_virilis",
                   "Drosophila_yakuba", "Fugu_rubripes", "Gallus_gallus",
                   "Homo_sapiens", "Monodelphis_domestica", "Mus_musculus",
                   "Pan_troglodytes", "Rattus_norvegicus", "SARS_coronavirus",
                   "Saccharomyces_cereviciae", "Saccharomyces_cerevisiae",
                   "Takifugu_rubripes", "Tetraodon_nigroviridis",
                   "Xenopus_tropicalis")
    wantedFiles <- c("refGene.txt.gz", "refLink.txt.gz", "cytoBand.txt.gz")
    accepts <-  paste(wantedFiles, collapse=",")
    wgetCmd <- paste("wget", "--timestamping", "--recursive", "--accept",
                     accepts)
    ## baseUrl must not end in a slash, as wget is sensitive to repeated /'s
    baseUrl <- gsub("/$", "", baseUrl)
    urls <- paste(baseUrl, organisms, "database", sep="/")
    wgetCmds <- paste(wgetCmd, urls)
    cmdRunner <- function(cmd) {
        cat("Running:", cmd, "\n")
        system(cmd)
    }
    returnCodes <- sapply(wgetCmds, cmdRunner)
    returnCodes
}


getPubDataGo <- function(baseUrl) {
    wget(baseUrl, levels=1, accepts="xml.gz")
}


getPubDataYeastGenome <- function(baseUrl) {
    ## XXX: trailing '/' matters for wget :-(
    litUrl <- paste(baseUrl, "literature_curation/", sep="/")
    chrUrl <- paste(baseUrl, "chromosomal_feature/", sep="/")
    wget(litUrl, levels=1)
    wget(chrUrl, levels=1)
}


getPubDataKegg <- function(baseUrl) {
    wget(baseUrl)
    cat("KEGG download complete, you'll need to untar the result\n")
}

getPubDataArabidopsis <- function(baseUrl) {
    TRUE
}
