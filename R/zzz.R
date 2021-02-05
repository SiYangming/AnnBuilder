.onLoad <- function(libname, pkgname) {
    require("methods", quietly=TRUE) || stop("Package methods unavailable!")
    require("utils", quietly=TRUE) || stop("Package utils unavailable!")
    require("Biobase", quietly = TRUE) || stop("Package Biobase unavailable!")
    require("XML", quietly = TRUE) || stop(paste("Package XML unavailable!",
                   "XML can be obtained from www.R-project.org CRAN"))
    require("annotate", quietly = TRUE) || stop("Package annoate unavailable")
    if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
            addVigs2WinMenu("AnnBuilder")
    }

    tmpEnv <- new.env(parent=NULL)
    data("pubDataURLs", "sourceURLs", package="AnnBuilder", envir=tmpEnv)
    options(AnnBuilderPublicDataUrls=tmpEnv$pubDataURLs)
    options(AnnBuilderSourceUrls=tmpEnv$sourceURLs)
    rm(tmpEnv)

}

