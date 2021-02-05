.onLoad <- function(libname, pkgname){
    options(show.error.messages=FALSE)
    getPfamDb()
    options(show.error.messages=TRUE)    
}

.onUnload <- function(libpath){
    closePfamDb()
}
