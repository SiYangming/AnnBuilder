## A master file for generating "autoloading" data package zzz.R files
## To setup for a new data package
## > library(Biobase)
## > copySubstitute("dataPkgLoad.R","<pkg>/R/zzz.R",list(PACKAGE=XXXXX))
## Then the file will be properly prepared

.getDatasets <- function(path) {
    rdas <- list.files(file.path(path, "data"),
                       pattern = "\\.*.rda")
    rdas <- gsub("(^.*)\.rda", "\\1"
                 , rdas)
    return(rdas)
}

.rdas <- .getDatasets(.find.package("@PACKAGE@"))
for (.i in .rdas)
    eval(substitute(obj <- delay((function(i, pkg){
        file <- file.path(system.file("data", package=pkg),
                          paste(i, ".rda", sep=""))
        env <- new.env()
        load(file, envir = env)
        get(i, env)
    })(obj, "@PACKAGE@")), list(obj=.i)))

rm(.getDatasets,.i,.rdas)

.no_S3_generics = TRUE
