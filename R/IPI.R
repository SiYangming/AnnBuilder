# Constructs an object of IPI that handles IPI data files
#

IPI <- function(srcUrl = getSrcUrl("IPI"),
                organism = organism,
                built,
                fromWeb = TRUE
                ){

    speciesName <- toupper(organism2species(organism))
    baseFile <- paste("ipi.", toupper(organism2species(organism)), ".dat.gz", sep="")
    
    new("IPI", srcUrl = srcUrl,
        baseFile = baseFile,
        built = ifelse(missing(built), getSrcBuilt("IPI"), built),
        fromWeb = fromWeb
        )
}

speciesNorganism <- function() {
    speciesNorganismTable <- rbind(
                                   c("human", "Homo sapiens"),
                                   c("mouse", "Mus musculus"),
                                   c("rat", "Rattus norvegicus"),
				   c("chick", "Gallus gallus"),
				   c("bovin", "Bos taurus")
                                   )
    colnames(speciesNorganismTable) <- c("species", "organism")
    return(speciesNorganismTable)
}

species2organism <- function(species) {
    speciesNorganismTable <- speciesNorganism()
    return(speciesNorganismTable[which(speciesNorganismTable[,1]==tolower(species)),2])
}

organism2species <- function(organism) {
    speciesNorganismTable <- speciesNorganism()
    return(speciesNorganismTable[which(speciesNorganismTable[,2]==organism),1])
}
