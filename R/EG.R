# Constructs an object of EG that handles Entrez Gene data files
#

EG <- function(srcUrl = getSrcUrl("EG"),
               parser = c(accession = file.path(.path.package("AnnBuilder"),
                      "scripts", "egAccParser")), baseFile = "",
               built, accession = "gene2accession.gz", info = "gene_info.gz",
               go = "gene2go.gz", pubmed = "gene2pubmed.gz",
               refseq = "gene2refseq.gz", unigene = "gene2unigene",
                            mim = "mim2gene", fromWeb = TRUE){
  
    new("EG", srcUrl = srcUrl, parser = parser, baseFile = baseFile,
        built = ifelse(missing(built), getSrcBuilt("EG"), built),
        accession = accession, info = info, go = go, pubmed = pubmed,
        refseq = refseq, unigene = unigene, mim = mim, fromWeb = fromWeb)
}



