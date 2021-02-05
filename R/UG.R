# Constructs an object of UG
#

UG <- function(srcUrl = getSrcUrl("UG", "Homo sapiens"),
               parser = file.path(.path.package("AnnBuilder"), "scripts",
                                                    "gbUGParser"),
               baseFile = "", organism = "human", built, fromWeb = TRUE){

    new("UG", srcUrl = srcUrl, parser = parser, baseFile = baseFile,
        organism = organism, built = ifelse(missing(built),
                             getSrcBuilt("UG", organism), built),
        fromWeb = fromWeb)
}



