HG <- function(srcUrl = " ftp://ftp.ncbi.nih.gov/pub/HomoloGene/old/hmlg.ftp",
                 fromWeb = TRUE, built = getSrcBuilt("hg")){

    new("HG", srcUrl = srcUrl, fromWeb = fromWeb, built = built)
}
