addNamespace <- function(pkgName, pkgPath,
                         hidePattern = c("QC", "MAPCOUNTS")){
##dataset should not be in a name space
#    fileNames <- gsub("\\.rda", "", list.files(file.path(pkgPath,
#                         pkgName, "data"), pattern = "rda"))
#    fileNames <- fileNames[-grep(paste(hidePattern, sep = "",
#                                      collapse = "|"), fileNames)]
#    write(paste("export(", paste(c(pkgName, fileNames), sep = "",
#                                 collapse = ", "), ")", sep = ""),
#          file = file.path(pkgPath, pkgName, "NAMESPACE"))
}

sealEnvs <- function(pkgName, pkgPath){
    fileNames <- list.files(file.path(pkgPath, pkgName, "data"),
                            pattern = "rda")
    for(i in fileNames){
        load(file.path(pkgPath, pkgName, "data", i))
        envName <- gsub("\\.rda", "", i)
        if(is.environment(get(envName))){
            tempEnv <- get(envName)
            lockEnvironment(tempEnv, binding = TRUE)
            assign(envName, tempEnv)
            save(list = envName, file = file.path(pkgPath, pkgName,
                                 "data", i))
        }
    }
}

