# Functions that return the correct URL that is going to be used to
# get annotation data from the desirable public data source.
#
# Copyright 2002, J. Zhang, all rights reserved.
#

.srcUrls <- function(name) {
    sourceURLs <- getOption("AnnBuilderSourceUrls")
    sourceURLs[[name]]
}


getSrcUrl <- function(src, organism = "Homo sapiens", xml = TRUE,
                      dateOnly = FALSE){
    switch(toupper(src),
           "LL" = return(getLLUrl()),
           "GP" = return(getUCSCUrl(organism)),
           "UG" = return(getUGUrl(organism)),
           "GO" = return(getGOUrl(xml, dateOnly)),
           "KEGG" = return(getKEGGUrl()),
           "GEO" = return(getGEOUrl()),
           "YG" = return(getYGUrl()),
           "HG" = return(getHGUrl()),
           "REFSEQ" = return(getRefSeqUrl(organism)),
           "EG" = return(getEGUrl()),
           "ALL" = return(getAllUrl(organism)),
           stop(paste("Source", src, "is not supported.",
                      "Has to be LL, GP, UG, GO, KEGG, ALL")))
}

getAllUrl <- function(organism){
    urls <- c(GP = getUCSCUrl(organism),
              UG = getUGUrl(organism), GO = getGOUrl(), KEGG = getKEGGUrl(),
              YG = getYGUrl(), HG = getHGUrl(), EG = getEGUrl())
    return(urls)
}


getEGUrl <- function(){
  return(.srcUrls("EG"))
}

getYGUrl <- function(){
    return(.srcUrls("YG"))
}

getLLUrl <- function(){
    return(.srcUrls("LL"))
}

getUCSCUrl <- function(organism, downloadSite){
    if (missing(downloadSite))
      downloadSite <- .srcUrls("GP")
    key <- switch(toupper(organism),
                  "HOMO SAPIENS" = "Homo_sapiens",
                  "MUS MUSCULUS" = "Mus_musculus",
                  "RATTUS NORVEGICUS" = "Rattus_norvegicus",
                  "DANIO RERIO" = "Danio_Rerio",
                  "CAENORHABDITIS ELEGANS" = "Caenorhabditis_elegans",
                  "DROSOPHILA MELANOGASTER" = "Drosophila_melanogaster",
                  NA)
    if(is.na(key)) {
        warning(paste("Organism", organism, "is not supported by GoldenPath (GP)."))
        return(NA)
    }
    gpUrl <- paste(downloadSite, key, "database/", sep="/")
    gpUrl
}


getUGUrl <- function(organism){
    #switch(toupper(organism),
    #       "HUMAN" = return(paste(.srcUrls("UG"), "/Hs.data.gz",
    #                        sep = "")),
    #       "MOUSE" = return(paste(.srcUrls("UG"), "/Mm.data.gz",
    #                        sep = "")),
    #       "RAT" = return(paste(.srcUrls("UG"),
    #                       "/Rn.data.gz", sep = "")),
    #       "ZEBRAFISH" = return(paste(.srcUrls("UG"), "/Dr.data.gz",
    #                      sep = "")),
    #       return(NA))
           #stop(paste("Organism", organism, "is not supported.")))

    return(paste(.srcUrls("UG"), "/", getUGShortName(organism),
                 ".data.gz", sep = ""))
}

# Figures out the url for the latest version of the GO data file
getGOUrl <- function(xml = TRUE, dateOnly = FALSE){
    ## XXX: We don't know the name of the file we want, so we
    ##      do screen scraping to extract the filename from
    ##      the HTML directory listing page
    con <- url(.srcUrls("GO"))
    dirIndexHtml <- readLines(con)
    close(con)
    goFilePat <- "href=\"go.*-termdb.*rdf-xml.gz\""
    goFile <- grep(goFilePat, dirIndexHtml, value=TRUE)
    goFilePat <- "^.*href=\"(go.*rdf-xml.gz)\".*$"
    goFile <- gsub(goFilePat, "\\1", goFile)
    paste(.srcUrls("GO"), goFile, sep="/")
}

getRefSeqUrl <- function(organism){
    switch(tolower(organism),
           "homo sapiens" = return("ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/"),
           "yeast" =
    return("ftp://ftp.ncbi.nih.gov/genomes/Saccharomyces_cerevisiae/"),
           stop(paste("Organism", organism, "not supported")))
}

getKEGGUrl <- function(){
    return(.srcUrls("KEGG"))
}

getGEOUrl <- function(){
    return(.srcUrls("GEO"))
}

getHGUrl <- function(){
    return(.srcUrls("HG"))
}

getUGShortName <- function(sciName){
    temp <- UGSciNames()
    return(names(temp[tolower(temp) == tolower(sciName)]))
}


UGSciNames <- function(){
    c(Aga = "Anopheles gambiae", Ame = "Apis mellifera",
      Bmo = "Bombyx mori", Bt = "Bos taurus",
      Cel = "Caenorhabditis elegans", Cfa = "Canis familiaris",
      Cin = "Ciona intestinalis", Gga = "Gallus allus",
      Dm = "Drosophila melanogaster", Dr = "Danio rerio",
      Hma = "Hydra magnipapillata",
      Oar = "Ovis aries", Ola = "Oryzias latipes",
      Omy = "Oncorhynchus mykiss",
      Sma = "Schistosoma mansoni", Ssa = "Salmo salar",
      Ssc = "Sus scrofa", Str = "Xenopus tropicalis",
      Xl = "Xenopus laevis", At = "Arabidopsis thaliana",
      Gma = "Glycine max", Han = "Helianthus annus",
      Hv = "Hordeum vulgare",  Lsa = " Lactuca sativa",
      Les = "Lycopersicon esculentum", Lco = "Lotus corniculatus",
      Mdo = "Malus domestica", Mtr = "Medicago truncatula",
      Os = "Oryza sativa", Pta = "Pinus taeda", Ptp = "Populus tremula",
      Ptp = "Populus tremuloides", Sbi = "Sorghum bicolor",
      Sof = "Saccharum officinarum", Ta = "Triticum aestivum",
      Vva = "Vitis vinifera", Zm = "Zea mays",
      Cre = "Chlamydomonas reinhardtii", Ddi = "Dictyostelium discoideum",
      Mgr = "Magnaporthe grisea", Ncr = "Neurospora crassa",
      Tgo = "Toxoplasma gondii", Hs = "Homo sapiens", Mm ="Mus musculus",
      Rn = "Rattus norvegicus")
}

getShortSciName <- function(sciName){
    temp <- unlist(strsplit(sciName, " +"))
    if(length(temp) < 2){
        warning(paste("Species name may not be a two-word scientic",
                      "name.", sciName, "was returned."))
    }
    return(paste(toupper(substr(temp[1], 1, 1)),
                 substr(temp[2], 1, 2), sep = ""))
}
