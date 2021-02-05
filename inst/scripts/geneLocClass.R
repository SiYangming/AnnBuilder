#### geneLoc objects are used by exonPkgBuilder to represent data
#### for the location of an exon
setClass("geneLoc", representation(probeId = "character",
                                chrom = "character",
                                strand = "character",
                                startLoc = "numeric",
                                endLoc = "numeric"))

# Set the get methods
if(!isGeneric("probeId")){
    setGeneric("probeId",
               function(object) standardGeneric("probeId"))
}
setMethod("probeId", "geneLoc",
          function(object) object@probeId)
if(!isGeneric("chrom")){
    setGeneric("chrom",
               function(object) standardGeneric("chrom"))
}
setMethod("chrom", "geneLoc",
          function(object) object@chrom)
if(!isGeneric("strand")){
    setGeneric("strand",
               function(object) standardGeneric("strand"))
}
setMethod("strand", "geneLoc",
          function(object) object@strand)
if(!isGeneric("startLoc")){
    setGeneric("startLoc",
               function(object) standardGeneric("startLoc"))
}
setMethod("startLoc", "geneLoc",
          function(object) object@startLoc)
if(!isGeneric("endLoc")){
    setGeneric("endLoc",
               function(object) standardGeneric("endLoc"))
}
setMethod("endLoc", "geneLoc",
          function(object) object@endLoc)
# Print the contents of an exon object
if(!isGeneric("print")){
    setGeneric("print",
               function(x, ...) standardGeneric("print"))
}
setMethod("print", "geneLoc",
          function(x, ...) {
               print("An object of class geneLoc")
               print(paste("Chromosome =", chrom(x)))
               print(paste("probe ID =", probeId(x)))
               print(paste("Strand =", strand(x)))
               print(paste("Beging Location =", startLoc(x)))
               print(paste("End location =", endLoc(x)))
          })

if(!isGeneric("inExon")){
    setGeneric("inExon",
               function(object, location, chr) standardGeneric("inExon"))
}
setMethod("inExon", "geneLoc",
          function(object, location, chr) {
              if(chr != chrom(object)){
                  return(NULL)
              }
              if(chr >= startLoc(object) && chr <= endLoc(object)){
                  return(probeId(object))
              }
              return(NULL)
          })
