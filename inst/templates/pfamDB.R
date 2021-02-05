pfam_globals <- new.env()

getPfamDb  <- function(){
    dbCon <- try(get("dbCon", pfam_globals))
    if(inherits(dbCon, "try-error") || is.null(dbCon)){
        SQLiteDriver <- dbDriver("SQLite")
        pfamSQLiteDbmsName <- system.file("data", "pfam.sqlite", package="PFAM")
        pfam_globals$dbCon <- dbConnect(SQLiteDriver, dbname=pfamSQLiteDbmsName)
        return(pfam_globals$dbCon)
    }else{
        return(dbCon)
    }
}

closePfamDb <- function(){
    dbDisconnect(pfam_globals$dbCon)
    pfam_globals$dbCon <- NULL
}

setPfamDb <- function(con){
    pfamCon <- try(get("dbCon", pfam_globals))
    if(!inherits(pfamCon, "try-error")){
        closePfamDb()
    }
    pfam_globals$dbCon <- con
}

getPfam <- function(pfamCon, FROM, TO, TABLE, FOR){
    if(length(FROM)!=1){
        stop("Argument \"FROM\" must be of length 1.")
    }

    TO2 <- paste(TO, sep="", collapse=", ")

    if(sum(duplicated(FOR))!=0){
        warning("Argument \"FOR\" is not unique!")
        FOR <- unique(FOR)
    }
    
    if(is.null(FOR)){
        query <- paste("select ", FROM, ", ", TO2, " from ", TABLE, ";", sep="")
    }else{
        FOR <- paste("\"", FOR, "\"", sep="", collapse=", ")
        query <- paste("select ", FROM, ", ", TO2, " from ", TABLE, " where ", FROM, " in (", FOR, ");", sep="")
    }
    tmpResult <- dbGetQuery(pfamCon, query)
    tmpResult[, FROM][is.na(tmpResult[, FROM])] <- "NA"
    result <- split(tmpResult[, TO], tmpResult[, FROM])
    if(length(TO)>1){
        result <- lapply(result, as.list)
    }
    return(result)
}
