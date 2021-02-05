.onLoad <- function(libname, pkgname) {
    packageStartupMessage(sprintf("\n*** Deprecation warning ***:\nThe package '%s' is deprecated and will not be supported after Bioconductor release 2.3.\nPlease switch to using AnnotationDbi for all of your annotation needs.\nIf you need to build an annotation package, please see the SQLForge vignette in the AnnotationDbi package.\n\n", pkgname))
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

    tmpEnv <- new.env()
    data("pubDataURLs", "sourceURLs", package="AnnBuilder", envir=tmpEnv)
    options(AnnBuilderPublicDataUrls=tmpEnv$pubDataURLs)
    options(AnnBuilderSourceUrls=tmpEnv$sourceURLs)
    rm(tmpEnv)
}
