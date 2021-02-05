# Constructors an object of pubRepo
#

pubRepo <- function(srcUrl, parser = "", baseFile = "",
                    built = "Not available", fromWeb = TRUE){
    new("pubRepo", srcUrl = srcUrl, parser = parser,
        baseFile = baseFile, built = built, fromWeb = TRUE)
}

