
PFAM <- function(srcUrl=getSrcUrl("PFAM"),
                 fromWeb=TRUE) {
    new("PFAM", srcUrl = srcUrl,
        ##built = ifelse(missing(built), getSrcBuilt("PFAM"), built),
        fromWeb = fromWeb
        )
}

pfamBuilder <- function(pkgName="PFAM", pkgPath, version, author,
                        fromWeb=TRUE, lazyLoad=TRUE, useTmp=FALSE, sqlFile=NULL){
    require("RSQLite", quietly=TRUE) || stop(paste("RSQLite is needed to build", pkgName, "package"))
    srcObjs <- list()
    srcObjs[["PFAM"]] <- PFAM(srcUrl=getSrcUrl("PFAM"), fromWeb=fromWeb)

    makeSrcInfo()
    createEmptyDPkg(pkgName, pkgPath, force = TRUE)
    cat(paste(pkgName, "package directory structure initialized\n"))

    if(is.null(sqlFile)){
        tableList <- parseData(srcObjs[["PFAM"]])
        pfamNamespace <- pfam2DB(tableList=tableList, pkgName=pkgName, pkgPath=pkgPath, useTmp=useTmp)
        cat("SQLite database file created\n")
    }else{
        cat("Copying", basename(sqlFile), "...")
        pfamNamespace <- NULL
        system(paste("cp", sqlFile, file.path(pkgPath, pkgName, "data", basename(sqlFile))))
        cat("Done\n")
        cat("SQLite database file created\n")
    }

    repList <- getRepList("PFAM", srcObjs)
    repList[["PKGNAME"]] <- pkgName

    cat("Creating documentations ...")
    ##writeDocs(baseName=NULL, pkgName, pkgPath, version, author,
    ##                  repList, "PKGNAME")

    writeDescription(pkgName, pkgPath, version, author)
    modifyDescription(pkgPath, pkgName)
    writePfamQC(tableList=tableList, pkgName=pkgName, pkgPath=pkgPath)
    
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates","GENERAL.Rd"),
                   file.path(pkgPath, pkgName, "man", paste(pkgName, ".Rd", sep = "")),
                   list(PKGNAME = pkgName, SOURCENBUILT =
                        paste("\n\t", repList, sep = "",
                              collapse = "\n")), "#")

    ##writeZZZ(pkgPath, pkgName)
    writeFun(pkgPath, pkgName)
    ##getDPStats(baseF=NULL, pkgName, pkgPath, isFile)
    copySubstitute(file.path(.path.package("AnnBuilder"),
                             "templates", "PKGNAMEQC.Rd"),
                   file.path(pkgPath, pkgName, "man",
                             paste(pkgName, "QC.Rd", sep = "")),
                   list(PKGNAME = pkgName), "#")
    ##copySubstitute(file.path(.path.package("AnnBuilder"),
    ##                         "templates", "PKGNAMEQCDATA.Rd"),
    ##               file.path(pkgPath, pkgName, "man",
    ##                         paste(pkgName, "QCDATA.Rd", sep = "")),
    ##               list(PKGNAME = pkgName), "#")
    cat(" Done\n")
    
    writeLines(readLines(file.path(.path.package("AnnBuilder"), "templates", "pfamZzz.R")),
               file.path(pkgPath, pkgName, "R", "zzz.R"))
    writeLines(readLines(file.path(.path.package("AnnBuilder"), "templates", "pfamDB.R")),
               file.path(pkgPath, pkgName, "R", "pfamDB.R"))
    ## won't export these functions on this stage
    ##pfamNamespace <- c(pfamNamespace, c("getPfamDb", "closePfamDb", "setPfamDb"))
    ##writeLines(readLines(file.path(.path.package("AnnBuilder"), "templates", "pfamDB.Rd")),
    ##           file.path(pkgPath, pkgName, "man", "pfamDB.Rd"))


    cat("Creating NAMESPACE file ...")
    ##pfamNamespace <- c(pfamNamespace, "getPfam")
    pfamNamespace <- c(pfamNamespace, "PFAM")
    pfamNamespace <- paste(pfamNamespace, sep="", collapse=", ")
    pfamNamespace <- paste("export(", pfamNamespace, ")", sep="")
    writeLines(pfamNamespace, file.path(pkgPath, pkgName, "NAMESPACE"))
    cat("Done\n")

    if (lazyLoad) {
        makeLLDB(file.path(pkgPath, pkgName), compress=TRUE)
    }
    
}

pfamParser <- function(pfamFile, fromWeb) {

    if(fromWeb) {
        tmpFile <- loadFromUrl(pfamFile)
    }else{
        unzipFile(fileName=basename(pfamFile),
                  where=dirname(pfamFile),
                  isgz = TRUE)
        tmpFile <- gsub("(.gz)", "\\", pfamFile)
    }
    keyDr <- getPfamDr(tmpFile)
    con <- file(tmpFile, "r")

    tableList <- list()
    ##for(key in c(keyDr, "basicInfo", "RM")){
    ##    tableList[[key]] <- NULL
    ##}
    
    while(TRUE) {
        ans <- readPfamRecord(con)
        ansAC <- getAC(ans$record)
        tableList[["basicInfo"]] <- rbind(tableList[["basicInfo"]],
                                          c(ansAC,
                                            getKeyValue(ans$record, key="#=GF ID   "),
                                            getKeyValue(ans$record, key="#=GF DE   "),
                                            getKeyValue(ans$record, key="#=GF TP   ")))
        
        ansRM <- getKeyValue(ans$record, key="#=GF RM   ")
        tableList[["RM"]] <- rbind(tableList[["RM"]], cbind(rep(ansAC, length(ansRM)), ansRM))

        for(key in setdiff(keyDr, c("SCOP", "PDB"))){
            keyValue <- getKeyValue(ans$record, key=paste("#=GF DR   ", key, sep=""), split="[ ;+]")
            tableList[[key]] <- rbind(tableList[[key]], cbind(rep(ansAC, length(keyValue)), keyValue))
        }
        
        keyValue <- getKeyValue(ans$record, key=paste("#=GF DR   ", "SCOP", sep=""), split="[ ;+]", which=1:2)
        tableList[["SCOP"]] <- rbind(tableList[["SCOP"]], cbind(rep(ansAC, nrow(keyValue)), keyValue))
        
        keyValue <- getKeyValue(ans$record, key=paste("#=GF DR   ", "PDB", sep=""), which=1:3)
        tableList[["PDB"]] <- rbind(tableList[["PDB"]], cbind(rep(ansAC, nrow(keyValue)), keyValue))
        
        if (ans$done) {
            break
        }
    }
    close(con)

    colnames(tableList[["basicInfo"]]) <- c("AC", "ID", "DE", "TP")
    tableList[["basicInfo"]] <- tableList[["basicInfo"]][apply(tableList[["basicInfo"]], 1, function(x) !all(is.na(x))), ]

    for(key in c("RM", setdiff(keyDr, c("SCOP", "PDB")))){
            colnames(tableList[[key]]) <- c("AC", key)
            tableList[[key]] <- tableList[[key]][apply(tableList[[key]], 1, function(x) !all(is.na(x))), ]

    }

    colnames(tableList[["SCOP"]]) <- c("AC", "SCOP", "PLACEMENT")
    tableList[["SCOP"]] <- tableList[["SCOP"]][apply(tableList[["SCOP"]], 1, function(x) !all(is.na(x))), ]

    colnames(tableList[["PDB"]]) <- c("AC", "PDB", "START", "END")
    tableList[["PDB"]] <- tableList[["PDB"]][apply(tableList[["PDB"]], 1, function(x) !all(is.na(x))), ]
    tableList[["PDB"]] <- matrix(gsub("^ ", "", tableList[["PDB"]]), nrow=nrow(tableList[["PDB"]]), ncol=ncol(tableList[["PDB"]]))

    return(tableList)
}

getPfamDr <- function(pfamFile){
    allDr <- system(paste("grep \"#=GF DR   \" ", pfamFile, sep=""), intern=TRUE)
    allDr <- gsub("#=GF DR   ", "", allDr)
    allDr <- sapply(allDr, function(x) unlist(strsplit(x, "[ ;+]"))[1])
    allDr <- unique(allDr)
    return(allDr)
}


readPfamRecord <- function(con) {
    buf <- NULL
    atEnd <- FALSE
    done <- FALSE
    while(!done) {
        line <- try(readLines(con, n=1, ok=FALSE), silent=TRUE)
        if (inherits(line, "try-error")) {
            atEnd <- TRUE
            done <- TRUE
        } else if (length(grep("^//$", line)) > 0) {
            done <- TRUE
        }
        buf <- c(buf, line)
    }
    list(record=buf, done=atEnd)
}

getKeyValue <- function(record, key, split=";", which=1) {
    keyValue <- grep(key, record, value=TRUE)
    if(length(keyValue)>0){
        keyValue <- gsub(key, "", keyValue)
        if(length(which)==1){
            keyValue <- sapply(strsplit(keyValue, split), function(x) x[x!=""][which])
        }else{
            keyValue <- lapply(strsplit(keyValue, split), function(x) x[x!=""][which])
            keyValue <- t(as.matrix(data.frame(keyValue)))
            colnames(keyValue) <- rownames(keyValue) <- NULL
        }
        return(keyValue)
    }else{
        if(length(which)==1){
            return(as.character(NA))
        }else{
            return(as.matrix(t(rep(as.character(NA), length(which)))))
        }
    }
}

getAC <- function(record){
    value <- getKeyValue(record, "#=GF AC")
    if(!is.na(value)) {
        value <- tail(unlist(strsplit(value, " +")), 1)
        value <- head(unlist(strsplit(value, "\\.")), 1)
        return(value)
    } else {
        return(as.character(NA))
    }
}

pfam2DB <- function(tableList, pkgName, pkgPath, useTmp=FALSE){
    require("RSQLite", quietly=TRUE) || stop(paste("RSQLite is needed to build", pkgName, "package."))

    pfamNamespace <- NULL
    
    SQLiteDriver <- dbDriver("SQLite")
    if(useTmp){
        pfamSQLiteDbmsName <- file.path(tempdir(), "pfam.sqlite")
    }else{
        pfamSQLiteDbmsName <- file.path(pkgPath, pkgName, "data", "pfam.sqlite")
    }
    pfamCon <- dbConnect(SQLiteDriver, dbname=pfamSQLiteDbmsName)

    writeMat2DB <- function(con, TABLE, mat) {
        mat <- matrix(gsub("'", "''", mat), nrow=nrow(mat), ncol=ncol(mat))
        mat <- data.frame(mat)
        sql <- NULL
        row2Sql <- function(row) {
            paste("insert into ", TABLE, " values (",
                  paste("'", row, "'", sep="", collapse=", "),
                  ");", sep="")
        }
        sql <- apply(mat, 1, row2Sql)
        sql <- paste(sql, collapse="\n ")
        dbSendQuery(con, sql)
    }
    
    cat("Creating basicInfo table ...")
    createTableSql <- "create table basicInfo (
                                AC varchar(7),  -- Accession number
                                ID varchar(15), -- Identification
                                DE varchar(80), -- Definition
                                TP varchar(6)   -- Type field
                                );"
    dbSendQuery(pfamCon, createTableSql)
    writeMat2DB(pfamCon, TABLE="basicInfo", mat=tableList$basicInfo)

    detailInfo <- getDetailInfo()

    for(tab in c("ID", "DE", "TP")){
        pfamNamespace <- writePfamMapFun(ONE="AC", TWO=tab, TABLE="basicInfo", pkgName=pkgName, pkgPath=pkgPath, pfamNamespace=pfamNamespace)
        writePfamMapDoc(ONE="AC", TWO=tab, TABLE="basicInfo", EXTRA="None", DETAILINFO=detailInfo[[tab]], pkgName=pkgName, pkgPath=pkgPath)
    }
    cat(" Done\n")
    
    for(tab in setdiff(names(tableList), c("basicInfo", "PDB", "SCOP"))){
        cat(paste("Creating ", tab, " table ...", sep=""))
        colLength <- apply(tableList[[tab]], 2, function(x) max(nchar(x)))
        createTableSql <- paste("create table ", tab, " ( AC varchar(7), ", tab, " varchar(", colLength[2], ") );", sep="")
        dbSendQuery(pfamCon, createTableSql)
        writeMat2DB(pfamCon, TABLE=tab, mat=tableList[[tab]])
        pfamNamespace <- writePfamMapFun(ONE="AC", TWO=tab, TABLE=tab, pkgName=pkgName, pkgPath=pkgPath, pfamNamespace=pfamNamespace)
        writePfamMapDoc(ONE="AC", TWO=tab, TABLE=tab, EXTRA="None", DETAILINFO=detailInfo[[tab]], pkgName=pkgName, pkgPath=pkgPath)
        cat(" Done\n")
    }

    cat(paste("Creating SCOP table ...", sep=""))
    colLength <- apply(tableList[["SCOP"]], 2, function(x) max(nchar(x)))
    createTableSql <- paste("create table ", "SCOP", " ( AC varchar(7), SCOP varchar(", colLength[2], "), PLACEMENT varchar(", colLength[3], ") );", sep="")
    dbSendQuery(pfamCon, createTableSql)
    writeMat2DB(pfamCon, TABLE="SCOP", mat=tableList[["SCOP"]])
    cat(" Done\n")
    
    cat(paste("Creating PDB table ...", sep=""))
    colLength <- apply(tableList[["PDB"]], 2, function(x) max(nchar(x)))
    createTableSql <- paste("create table ", "PDB", " ( AC varchar(7), PDB varchar(", colLength[2], "), startPoint int, endPoint int );", sep="")
    dbSendQuery(pfamCon, createTableSql)
    writeMat2DB(pfamCon, TABLE="PDB", mat=tableList[["PDB"]])
    cat(" Done\n")

    ## a new way to copy file :-)
    writeLines(readLines(file.path(.path.package("AnnBuilder"), "templates", "pfamSpeMap.R")),
               file.path(pkgPath, pkgName, "R", "pfamSpeMap.R"))
    pfamNamespace <- c(pfamNamespace, c("pfamAC2SCOP", "pfamSCOP2AC", "pfamAC2PDB", "pfamPDB2AC"))
    writePfamMapDoc(ONE="AC", TWO="PDB", TABLE="PDB", EXTRA="Start and end points", DETAILINFO=detailInfo[["PDB"]], pkgName=pkgName, pkgPath=pkgPath)
    writePfamMapDoc(ONE="AC", TWO="SCOP", TABLE="SCOP", EXTRA="Placement", DETAILINFO=detailInfo[["SCOP"]], pkgName=pkgName, pkgPath=pkgPath)

    if(useTmp){
        cat("Copying SQLite file ...")
        if(.Platform$OS.type == "unix"){
            system(paste("cp ", pfamSQLiteDbmsName, " ", pkgPath, "/", pkgName, "/data/", sep=""))
            cat(" Done\n")
        }else if(.Platform$OS.type == "windows"){
            system(paste("copy ", pfamSQLiteDbmsName, " ", pkgPath, "/", pkgName, "/data/", sep=""))
            cat(" Done\n")                        
        }else{
            cat(" ERROR\n")
            stop("Unable to create RSQLite database file.  Please set the argument useTmp=FALSE .")
        }
    }

    ## error occurred when trying to disconnect pfamCon.  Remove it
    ##dbDisconnect(pfamCon)
    cat("All tables are created.\n")
    return(pfamNamespace)
}

writePfamMapFun <- function(ONE, TWO, TABLE, pkgName, pkgPath, pfamNamespace){
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates", "pfam.temp"),
                   file.path(pkgPath, pkgName, "R", paste("pfam", gsub("_", "", ONE), "2", gsub("_", "", TWO), ".R", sep="")),
                   list(FROM=ONE,
                        TO=TWO,
                        FROMNAME=(gsub("_", "", ONE)),
                        TONAME=(gsub("_", "", TWO)),
                        FOR=tolower(gsub("_", "", ONE)),
                        TABLE=TABLE
                        ),
                   "#")
    pfamNamespace <- c(pfamNamespace, paste("pfam", gsub("_", "", ONE), "2", gsub("_", "", TWO), sep=""))
    
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates", "pfam.temp"),
                   file.path(pkgPath, pkgName, "R", paste("pfam", gsub("_", "", TWO), "2", gsub("_", "", ONE), ".R", sep="")),
                   list(FROM=TWO,
                        TO=ONE,
                        FROMNAME=(gsub("_", "", TWO)),
                        TONAME=(gsub("_", "", ONE)),
                        FOR=tolower(gsub("_", "", TWO)),
                        TABLE=TABLE
                        ),
                   "#")
    pfamNamespace <- c(pfamNamespace, paste("pfam", gsub("_", "", TWO), "2", gsub("_", "", ONE), sep=""))
    return(pfamNamespace)
}

writePfamMapDoc <- function(ONE, TWO, TABLE, EXTRA, DETAILINFO, pkgName, pkgPath){

    if(is.null(DETAILINFO)){
        DETAILINFO <- "None"
    }
    
    tmpRead <- readLines("ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/userman.txt")
    
    longName <- function(shortName){
        whatWeWant <- grep(shortName, tmpRead, value=TRUE)
        if(length(whatWeWant)==0){
            return(paste(shortName, "ID"))
        }else{
            whatWeWant <- unlist(strsplit(whatWeWant, ":"))[1]
            whatWeWant <- unlist(strsplit(whatWeWant, "   "))[2]
            if(whatWeWant==""){
                return(paste(gsub("_", "\\\\_", shortName, fixed=TRUE), "ID"))
            }else{
                return(whatWeWant)
            }
        }
    }
    
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates", "pfam.Rd"),
                   file.path(pkgPath, pkgName, "man", paste("pfam", gsub("_", "", ONE), "2", gsub("_", "", TWO), ".Rd", sep="")),
                   list(FROM=gsub("_", "", ONE),
                        TO=gsub("_", "", TWO),
                        FOR=tolower(gsub("_", "", ONE)),
                        TABLE=TABLE,
                        EXTRA=EXTRA,
                        DETAILINFO=DETAILINFO,
                        FROMNAME=longName(ONE),
                        TONAME=longName(TWO)
                        ),
                   "#")
    
    copySubstitute(file.path(.path.package("AnnBuilder"), "templates", "pfam.Rd"),
                   file.path(pkgPath, pkgName, "man", paste("pfam", gsub("_", "", TWO), "2", gsub("_", "", ONE), ".Rd", sep="")),
                   list(FROM=gsub("_", "", TWO),
                        TO=gsub("_", "", ONE),
                        FOR=tolower(gsub("_", "", TWO)),
                        TABLE=TABLE,
                        EXTRA="None",
                        DETAILINFO=DETAILINFO,
                        FROMNAME=longName(TWO),
                        TONAME=longName(ONE)
                        ),
                   "#")
}

modifyDescription <- function(pkgPath, pkgName){
    tmpRead <- readLines(file.path(pkgPath, pkgName, "DESCRIPTION"))
    whatWeWant <- grep("^Depends:", tmpRead)
    tmpRead[whatWeWant] <- paste(tmpRead[whatWeWant], ", RSQLite", sep="")
    writeLines(tmpRead, file.path(pkgPath, pkgName, "DESCRIPTION"))
}

writePfamQC <- function(tableList, pkgName, pkgPath){

    QC <- paste("\n\nQuality control information for ", pkgName, "\n", sep="")
    QC <- paste(QC, "Date built: ", getDate(pkgName, pkgPath), "\n", sep="")
    
    numMapped <- sapply(tableList,
                        function(x){
                            tmpSplit <- split(x[, 2], x[, 1])
                            return(sum(sapply(tmpSplit, function(y) !all(is.na(y)))))
                        }
                        )
    QC <- paste(QC, "Number of Pfam entry: ", numMapped[1], "\n", sep="")
    QC <- paste(QC, "Mappings found from Pfam accession numbers: ", "\n", sep="")

    ACbasicInfo <- setdiff(colnames(tableList[["basicInfo"]]), "AC")
    mappings <- sort(c(setdiff(names(numMapped), "basicInfo"), ACbasicInfo))
    
    for(map in mappings){
        if(map %in% ACbasicInfo){
            numMap <- numMapped[1]
        }else{
            numMap <- numMapped[map]
        }
        QC <- paste(QC, "\tpfamAC2", gsub("_", "", map), " found ", numMap, " of ", numMapped[1], "\n", sep="")
    }
    
    numRevMapped <- sapply(tableList,
                           function(x){
                               tmpSplit <- split(x[, 1], x[, 2])
                               return(sum(sapply(tmpSplit, function(y) !all(is.na(y)))))
                           }
                           )

    numRevMapped <- c(numRevMapped, sapply(ACbasicInfo,
                                           function(x){
                                               tmpSplit <- split(tableList[["basicInfo"]][, 1], tableList[["basicInfo"]][, x])
                                               return(sum(sapply(tmpSplit, function(y) !all(is.na(y)))))
                                           }
                                           )
                      )
    
    QC <- paste(QC, "Reverse mappings found to Pfam accession numbers: ", "\n", sep="")
    for(map in setdiff(sort(names(numRevMapped)), "basicInfo")){
        QC <- paste(QC, "\tpfam", gsub("_", "", map), "2AC found ", numRevMapped[map], "\n", sep="")
    }

    dataPath <- file.path(pkgPath, pkgName, "data")
    qcVar <- paste(pkgName, "QC", sep = "")
    qcFile <- file.path(dataPath, paste(qcVar, "rda", sep="."))
    assign(qcVar, QC)
    save(list=qcVar, file=qcFile)
}

getDetailInfo <- function(){
    detailInfo <- list()
    detailInfo[["CAZY"]] <- "The CAZy database (\\\\url{http://afmb.cnrs-mrs.fr/CAZY/}) describes the families of structurally-related catalytic and carbohydrate-binding modules (or functional domains) of enzymes that degrade, modify, or create glycosidic bonds."
    detailInfo[["HOMSTRAD"]] <- "HOMSTRAD (HOMologous STRucture Alignment Database, \\\\url{http://www-cryst.bioc.cam.ac.uk/homstrad/}) is a curated database of structure-based alignments for homologous protein families.  Reference: Mizuguchi K, Deane CM, Blundell TL, Overington JP. (1998) HOMSTRAD: a database of protein structure alignments for homologous families. Protein Science 7:2469-2471."
    detailInfo[["INTERPRO"]] <- "\\\\url{http://www.ebi.ac.uk/interpro/}"
    detailInfo[["LOAD"]] <- "None"
    detailInfo[["MEROPS"]] <- "The MEROPS database (\\\\url{http://merops.sanger.ac.uk/}) is an information resource for peptidases (also termed proteases, proteinases and proteolytic enzymes) and the proteins that inhibit them.  Reference: Rawlings, N.D., Tolle, D.P. & Barrett, A.J. (2004) MEROPS: the peptidase database. Nucleic Acids Res. 32 Database issue, D160-D164"
    detailInfo[["MIM"]] <- "MIM (a.k.a. OMIM, \\\\url{http://www.ncbi.nlm.nih.gov/omim/}) is a catalog of human genes and genetic disorders authored and edited by Dr. Victor A. McKusick and his colleagues at Johns Hopkins and elsewhere.  Reference:  <MIM> MIM: McKusick, V.A.: Mendelian Inheritance in Man. A Catalog of Human Genes and Genetic Disorders. Baltimore: Johns Hopkins University Press, 1998 (12th edition). <OMIM> Online Mendelian Inheritance in Man, OMIM (TM). McKusick-Nathans Institute for Genetic Medicine, Johns Hopkins University (Baltimore, MD) and National Center for Biotechnology Information, National Library of Medicine (Bethesda, MD), 2000"
    detailInfo[["PDB"]] <- "PDB (\\\\url{http://www.rcsb.org/pdb/index.html}), the single worldwide repository for the processing and distribution of 3-D biological macromolecular structure data.  Reference: H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne:  The Protein Data Bank.   Nucleic Acids Research  ,  28  pp. 235-242 (2000)"
    detailInfo[["PFAMB"]] <- "Pfam-B is an automatic clustering of the rest of SWISSPROT and TrEMBL derived from the PRODOM database."
    detailInfo[["PRINTS"]] <- "PRINTS (\\\\url{http://umber.sbs.man.ac.uk/dbbrowser/PRINTS/}) is a compendium of protein fingerprints."
    detailInfo[["PROSITE"]] <- "PROSITE (\\\\url{http://www.expasy.ch/prosite/}) is a database of protein families and domains. It consists of biologically significant sites, patterns and profiles that help to reliably identify to which known protein family (if any) a new sequence belongs.  Reference: Hulo N., Sigrist C.J.A., Le Saux V., Langendijk-Genevaux P.S., Bordoli L., Gattiker A., De Castro E., Bucher P., Bairoch A. Recent improvements to the PROSITE database. Nucl. Acids. Res. 32:D134-D137(2004)"
    detailInfo[["RM"]] <- "Reference Medline (\\\\url{http://www.ncbi.nlm.nih.gov/PubMed/})"
    detailInfo[["SCOP"]] <- "Structural Classification of Proteins (\\\\url{http://scop.mrc-lmb.cam.ac.uk/scop/index.html}).  Reference: Murzin A. G., Brenner S. E., Hubbard T., Chothia C. (1995). SCOP: a structural classification of proteins database for the investigation of sequences and structures. J. Mol. Biol. 247, 536-540"
    detailInfo[["SMART"]] <- "SMART (a Simple Modular Architecture Research Tool, \\\\url{http://smart.embl-heidelberg.de/}) allows the identification and annotation of genetically mobile domains and the analysis of domain architectures.  Reference: (1) Schultz et al. (1998) Proc. Natl. Acad. Sci. USA 95, 5857-5864.  (2) Letunic et al. (2004) Nucleic Acids Res 32, D142-D144"
    detailInfo[["TC"]] <- "None"
    detailInfo[["URL"]] <- "None"    
    return(detailInfo)
}
