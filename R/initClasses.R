# A generic class that reads or downloads data from public data repositories.
# srcUrl -  the url for the cgi script that initiates a query against
#            the databse. The value at the time of coding is:
#            "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
#
    setClass("pubRepo", representation(srcUrl = "character",
                                       parser = "character",
                                       baseFile = "character",
                                       built = "character",
                                       fromWeb = "logical"))

    # Set the get methods
    if(!isGeneric("srcUrl")){
        setGeneric("srcUrl",
                   function(object) standardGeneric("srcUrl"))
    }
    setMethod("srcUrl", "pubRepo",
              function(object) object@srcUrl)

    if(!isGeneric("builtInfo")){
        setGeneric("builtInfo",
                   function(object) standardGeneric("builtInfo"))
    }
    setMethod("builtInfo", "pubRepo",
              function(object) object@built)
    if(!isGeneric("fromWeb")){
        setGeneric("fromWeb",
                   function(object) standardGeneric("fromWeb"))
    }
    setMethod("fromWeb", "pubRepo",
              function(object) object@fromWeb)
    if(!isGeneric("fromWeb<-")){
        setGeneric("fromWeb<-", function(object, value)
                   standardGeneric("fromWeb<-"))
    }
    setReplaceMethod("fromWeb", "pubRepo", function(object, value){
                  object@fromWeb <- value; object})
    # Define the replace methods
    if(!isGeneric("srcUrl<-")){
        setGeneric("srcUrl<-", function(object, value)
                   standardGeneric("srcUrl<-"))
    }
    setReplaceMethod("srcUrl", "pubRepo", function(object, value){
                  object@srcUrl <- value; object})
    if(!isGeneric("baseFile")){
        setGeneric("baseFile",
                   function(object) standardGeneric("baseFile"))
    }
    setMethod("baseFile", "pubRepo",
              function(object) object@baseFile)

    # Define the replace methods
    if(!isGeneric("baseFile<-")){
        setGeneric("baseFile<-", function(object, value)
                   standardGeneric("baseFile<-"))
    }
    setReplaceMethod("baseFile", "pubRepo", function(object, value){
                  object@baseFile <- value; object})
    if(!isGeneric("parser")){
        setGeneric("parser",
                   function(object) standardGeneric("parser"))
    }
    setMethod("parser", "pubRepo",
              function(object) object@parser)
    if(!isGeneric("parser<-")){
        setGeneric("parser<-", function(object, value)
                   standardGeneric("parser<-"))
    }
    setReplaceMethod("parser", "pubRepo", function(object, value){
                  object@parser <- value; object})
    # Defines functions
    if(!isGeneric("readData")){
        setGeneric("readData",
                   function(object, ...)
                   standardGeneric("readData"))
    }
    setMethod("readData", "pubRepo",
              function(object, ...){
                  if(fromWeb(object)){
                      conn <- url(srcUrl(object))
                  }else{
                      conn <- file(srcUrl(object))
                  }
                  temp <- readLines(conn)
                  close(conn)
                  return(temp)})
    if(!isGeneric("downloadData")){
        setGeneric("downloadData",
                   function(object, dist)
                   standardGeneric("downloadData"))
    }
    setMethod("downloadData", "pubRepo",
              function(object, dist)
                  return(loadFromUrl(srcUrl(object), dist)))
    if(!isGeneric("parseData")){
        setGeneric("parseData", function(object, ...)
                   standardGeneric("parseData"))
    }
    setMethod("parseData", "pubRepo", function(object, sep = "\t",
                                               ncol = 2, mergeKey = TRUE){
        if(fromWeb(object)){
            srcData <- downloadData(object, "")
        }else{
            srcData <- srcUrl(object)
        }
        tempOut <- tempfile("tempOut")
        obtained <- matrix(scan(fileMuncher(tempOut, baseFile(object),
                           srcData, parser(object)), what = "character",
                           sep = sep, quote = "", quiet = TRUE,
                           strip.white = TRUE, comment.char = ""),
                           ncol = ncol, byrow = TRUE)
        if(fromWeb(object)){
            unlink(srcData)
        }
        unlink(tempOut)
        if(nrow(obtained) <= 1 || !mergeKey){
            return(obtained)
        }else{
            return(mergeRowByKey(obtained))
        }})


# Sub class of pubRepo that reads/downloads data from GEO
# srcUrl -  the url for the cgi script that initiates a query against
#            the databse. The value at the time of coding is:
#            "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
    setClass("GEO", contains = "pubRepo")
    # Redefines readData
    setMethod("readData", "GEO",
              function(object, GEOAccNum)
              return(queryGEO(object, GEOAccNum)))


# Sub class of pubRepo that reads/downloads data from yeast genomic
# data.
# srcUrl - the base url for the ftp site for downloading. The value at
#          the time of coding is:
#          "ftp://genome-ftp.stanford.edu/pub/yeast/data_download/"
     setClass("YG", contains = "pubRepo")
     # Redefines readData
     setMethod("readData", "YG",
              function(object, extenName, cols2Keep, sep)
              return(getYeastData(srcUrl(object), extenName,
                                      cols2Keep, sep)))
# Sub class of pubRepo that processes GO data

    setClass("GO", contains = "pubRepo")
    setMethod("readData", "GO",
              function(object, xml = TRUE){
                  if(fromWeb(object)){
                      fileName <- loadFromUrl(srcUrl(object), "")
                  }else{
                      fileName <- srcUrl(object)
                  }
                  if(xml){
                      goData <- GOXMLParser(fileName)
                  }else{
                      goData <- readLines(fileName)
                  }
                  unlink(fileName)
                  return(goData)})

    setClass("LL", contains = "pubRepo")

    setClass("EG", contains = "pubRepo",
             representation(accession = "character",
                            info = "character",
                            go = "character",
                            pubmed = "character",
                            refseq = "character",
                            unigene = "character",
                            mim = "character"))
    setMethod("builtInfo", "pubRepo",
              function(object) object@built)

    setMethod("parseData", "EG", function(object, what = "gene2accession.gz",
                                    sep = "\t", ncol = 2, mergeKey = TRUE){
        if(fromWeb(object)){
          srcData <- loadFromUrl(paste(srcUrl(object), what, sep = "/"))
        }else{
          srcData <- srcUrl(object)
        }
        tempOut <- tempfile("tempOut")
        obtained <- matrix(scan(fileMuncher(tempOut, baseFile(object),
                                srcData, parser(object)), what = "character",
                                sep = sep, quote = "", quiet = TRUE,
                                strip.white = TRUE, comment.char = ""),
                           ncol = ncol, byrow = TRUE)
        if(fromWeb(object)){
          unlink(srcData)
        }
        unlink(tempOut)
        if(nrow(obtained) <= 1 || !mergeKey){
          return(obtained)
        }else{
          return(mergeRowByKey(obtained))
        }
    })

    setClass("UG", contains = "pubRepo",
             representation(organism = "character"))
    if(!isGeneric("orgName")){
        setGeneric("orgName",
                   function(object) standardGeneric("orgName"))
    }
    setMethod("orgName", "UG",
              function(object) object@organism)
    if(!isGeneric("orgName<-")){
        setGeneric("orgName<-", function(object, value)
                   standardGeneric("orgName<-"))
    }
    setReplaceMethod("orgName", "UG", function(object, value){
                  object@organism <- value; object})

    setClass("KEGG", contains = "UG")
    # Defines specific functions
    if(!isGeneric("findIDNPath")){
        setGeneric("findIDNPath", function(object)
                   standardGeneric("findIDNPath"))
    }
    # Get a named vector that labels KEGG pathway names with pathway ids.
    setMethod("findIDNPath", "KEGG", function(object)
              return(getKEGGIDNName(object)))
    if(!isGeneric("mapLL2ECNPName")){
        setGeneric("mapLL2ECNPName", function(object)
                   standardGeneric("mapLL2ECNPName"))
    }
    setMethod("mapLL2ECNPName", "KEGG", function(object)
              return(getLLPathMap(srcUrl(object), findIDNPath(object),
                               orgName(object), fromWeb = fromWeb(object))))

    setClass("GP", contains = "UG")
    # Defines specific functions
    if(!isGeneric("getStrand")){
        setGeneric("getStrand", function(object, ...)
                   standardGeneric("getStrand"))
    }
    # Parse refLink and refGene data files to get chromosomal location data
    setMethod("getStrand", "GP", function(object)
              return(getChroLocation(srcUrl(object),
                                     fromWeb = fromWeb(object))))

# Sub class of pubRepo that reads/downloads data from Homologene
# srcUrl -  the url for the file containing homology data
# ftp://ftp.ncbi.nih.gov/pub/HomoloGene/hmlg.ftp
    setClass("HG", contains = "pubRepo")
    # Redefines readData
    setMethod("readData", "HG",
              function(object)
              return(procHomoData(srcUrl(object))))
