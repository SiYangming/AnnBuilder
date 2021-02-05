makeLLDB <- function(packageDir, compress=TRUE) {
    dataDir <- file.path(packageDir, "data")
    dataEnv <- new.env(hash=TRUE)
    tmpEnv <- new.env()
    f0 <- files <- tools:::list_files_with_type(dataDir, "data")
    files <- unique(basename(files))
    dlist <- vector("list", length(files))
    names(dlist) <- files
    loaded <- character(0)
    for(i in 1:length(f0)) {
        load(file = f0[i], envir = dataEnv)
        load(file = f0[i], envir = tmpEnv)
        tmp <- ls(envir = tmpEnv, all.names = TRUE)
        rm(list = tmp, envir = tmpEnv)
        dlist[[files[i]]] <- tmp
        loaded <- c(loaded, tmp)
    }
    dup<- duplicated(loaded)
    if(any(dup))
      warning("object(s) ", paste(sQuote(loaded[dup]), collapse=", "),
              " are created by more than one data call")

    if(length(loaded)) {
        dbbase <- file.path(dataDir, "Rdata")
        tools:::makeLazyLoadDB(dataEnv, dbbase, compress = compress)
        .saveRDS(dlist, file.path(dataDir, "Rdata.rds"))
        cat("The following data sets have been added to the database and will be removed:\n")
        print(f0)
        unlink(f0)
        if (any(file.exists(f0)))
          warning("Some rda files were not removed, be sure to check the data dir before building")
        if(file.exists(file.path(dataDir, "filelist")))
          unlink(file.path(dataDir, c("filelist", "Rdata.zip")))
    }
}
