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
    mirrorUrls[idx] <- "kegg/pathway/organisms"
    df <- data.frame(name=names(pubUrls), url=mirrorUrls)
    write.table(file=file, df, row.names=FALSE, quote=FALSE)
}


readSourceUrlConfig <- function(file, urlPrefix) {
    df <- read.table(file, colClasses=rep("character", 2), header=TRUE, stringsAsFactors=FALSE)
    urls = df$url
    urls <- as.list(urls)
    if (!missing(urlPrefix))
      urls <- paste(urlPrefix, urls, sep="/")
    urls <- as.list(urls)
    names(urls) = df$name
    urls
}

    
downloadSourceData <- function(passive=FALSE) {
    pubDataURLs <- getOption("AnnBuilderPublicDataUrls")
    getPubDataHomoloGene(pubDataURLs[["HG"]], passive)
    getPubDataLocusLink(pubDataURLs[["LL"]], passive)
    getPubDataUniGene(pubDataURLs[["UG"]], passive)
    getPubDataEntrezGene(pubDataURLs[["EG"]], passive)
    getPubDataGoldenPath(pubDataURLs[["GP"]], passive)
    getPubDataGo(pubDataURLs[["GO"]], passive)
    getPubDataYeastGenome(pubDataURLs[["YG"]], passive)
    getPubDataKegg(pubDataURLs[["KEGG"]], passive)
    getPubDataArabidopsis(pubDataURLs[["AT"]], passive)
    getPubDataIpi(pubDataURLs[["IPI"]], passive)
    getPubDataCMAP(pubDataURLs[["CMAP"]], passive)
    getPubDataYEAST(pubDataURLs[["YEAST"]], passive)
    getPubDataKEGGGENOME(pubDataURLs[["KEGGGENOME"]], passive)
    getPubDataPfam(pubDataURLs[["PFAM"]], passive)

    if(.Platform$OS.type == "unix") {
        ## Unpack KEGG
        dir.create("kegg/pathway/organisms", recursive=TRUE)
        oldWD <- getwd()
        ## Remove 'http://' to get downloaded tarball name
        keggUrl <- pubDataURLs[["KEGG"]]
        keggTarball <- substr(keggUrl, 7, nchar(keggUrl))
        keggTarball <- file.path(oldWD, keggTarball)

        setwd(file.path(oldWD, "kegg", "pathway", "organisms"))
        system(paste("tar", "xzf", keggTarball))
	system("mv map_title.tab ../.")
        setwd(oldWD)
    } else {
        cat("Don't forget to untar KEGG!\n")
    }

}


getWgetCmd <- function(url, levels, accepts, passive=FALSE) {
    if (!missing(accepts)) {
        accepts <- paste(accepts, collapse=",")
        accepts <- paste("--accept", accepts)
    } else {
        accepts <- ""
    }
    if (missing(levels))
      levels <- 5
    levels <- paste("--level", levels)
    if (passive)
      passive <- "--passive-ftp"
    else
      passive <- ""
    
    globalOpts <- paste("--recursive", "--timestamping", "--retr-symlinks", 
                          levels, accepts, passive)

    wgetCmd <- paste("wget", globalOpts, url)
    wgetCmd
}


wget <- function(url, levels, accepts, passive=FALSE) {
    wgetCmd <- getWgetCmd(url, levels, accepts, passive)
    cat("Running:", wgetCmd, "\n")
    system(wgetCmd)
}


getPubDataHomoloGene <- function(baseUrl, passive) {
    wget(baseUrl, levels=1, accepts="hmlg.ftp", passive)
    currentUrl <- gsub("old", "current", baseUrl)
    wget(currentUrl, levels=1, accepts="hmlg.ftp", passive)
}


getPubDataLocusLink <- function(baseUrl, passive) {
    wget(baseUrl, levels=2, accepts="LL_tmpl.gz", passive)
}


getPubDataUniGene <- function(baseUrl, passive) {
    wget(baseUrl, levels=2, accepts=c(".data.gz", ".info"), passive)
}


getPubDataEntrezGene <- function(baseUrl, passive) {
    wget(baseUrl, levels=1, passive=passive)
}


getPubDataGoldenPath <- function(baseUrl, passive) {
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
    ## baseUrl must not end in a slash, as wget is sensitive to repeated /'s
    baseUrl <- gsub("/$", "", baseUrl)
    urls <- paste(baseUrl, organisms, "database", sep="/")
    wgetCmds <- getWgetCmd(url=urls, accepts=accepts, passive=passive)
    ##wgetCmds <- paste(wgetCmd, urls)
    cmdRunner <- function(cmd) {
        cat("Running:", cmd, "\n")
        system(cmd)
    }
    returnCodes <- sapply(wgetCmds, cmdRunner)
    returnCodes
}


getPubDataGo <- function(baseUrl, passive) {
    wget(baseUrl, levels=1, accepts="xml.gz", passive)
}


getPubDataYeastGenome <- function(baseUrl, passive) {
    ## XXX: trailing '/' matters for wget :-(
    litUrl   <- paste(baseUrl, "literature_curation/", sep="/")
    chrUrl   <- paste(baseUrl, "chromosomal_feature/", sep="/")
    regisUrl <- paste(baseUrl, "gene_registry/",       sep="/")
    wget(litUrl,   levels=1, passive=passive)
    wget(chrUrl,   levels=1, passive=passive)
    wget(regisUrl, levels=1, passive=passive)
}


getPubDataKegg <- function(baseUrl, passive) {
    wget(baseUrl, passive=passive)
    cat("KEGG download complete, you'll need to untar the result\n")
}

getPubDataArabidopsis <- function(baseUrl, passive) {
    LocusUrl <- paste(baseUrl, "Microarrays/Affymetrix/", sep="/")
    GenesUrl <- paste(baseUrl, "Genes/", sep="/")
    EstUrl <- paste(baseUrl, "Genes/est_mapping/", sep="/")
    GOUrl <- paste(baseUrl, "Ontologies/Gene_Ontology/", sep="/")
    PathwaysUrl <- paste(baseUrl, "Pathways/", sep="/")
    POUrl <- paste(baseUrl, "User_Requests/", sep="/")
    wget(LocusUrl, levels=1, passive=passive)
    wget(GenesUrl, levels=1, passive=passive)
    wget(EstUrl, levels=1, passive=passive)
    wget(GOUrl, levels=1, passive=passive)
    wget(PathwaysUrl, levels=1, passive=passive)
    wget(POUrl, levels=1, passive=passive)
}

getPubDataIpi <- function(baseUrl, passive) {
    wget(baseUrl, levels=1, accepts="dat.gz", passive)
}

getPubDataCMAP <- function(baseUrl, passive) {
    wget(baseUrl, levels=1, accepts="xml.gz", passive)
}

getPubDataYEAST <- function(baseUrl, passive) {
    wget(baseUrl, levels=1, passive=passive, accepts="tab")
}

getPubDataKEGGGENOME <- function(baseUrl, passive) {
    wget(baseUrl, levels=1, passive=passive)
}

getPubDataPfam <- function(baseUrl, passive){
    wget(baseUrl, passive=passive)
}
