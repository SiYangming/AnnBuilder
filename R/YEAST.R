
YEAST <- function(srcUrl=getSrcUrl("YEAST"),
                  baseFile="domains.tab",
                  built,
                  fromWeb=TRUE) {
    new("YEAST", srcUrl = srcUrl,
        baseFile = baseFile,
        built = ifelse(missing(built), getSrcBuilt("YEAST"), built),
        fromWeb = fromWeb
        )
}
