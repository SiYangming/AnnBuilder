# This function writes header information to an XML file.
#
# Copyright 2002 J. Zhang, all rights reserved.
#
writeXMLHeader <- function (outName, fileCol, name, version,
                            organism = "human"){

    writeLine <- function(toWrite){
        write(x = toWrite, file = outName, append = TRUE)
    }

    writeLine(paste("<?xml version = \"1.0\" encoding = ",
                       "\"UTF-8\" standalone = \"yes\"?>", sep = ""))

    writeLine(paste("<!DOCTYPE AnnBuilder: SYSTEM ",
                        "\"http://www.bioconductor.org/datafiles/dtds/",
                        "annotate.dtd\">", sep = ""))

    writeLine(paste("<AnnBuilder:Annotate xmlns:AnnBuilder = ",
                       "'http://www.bioconductor.org/AnnBuilder/'>",
                       sep = ""))

    writeLine("<AnnBuilder:Attr>")
    writeLine(paste("<AnnBuilder:Target value = \"", name,
                    "\"/>", sep = ""))
    writeLine(paste("<AnnBuilder:DateMade value = \"", date(),
                    "\"/>", sep = ""))
    writeLine(paste("<AnnBuilder:Version value = \"", version,
                    "\"/>", sep = ""))

    for (i in getDSrc(organism)){
        statement <- paste("<AnnBuilder:SourceFile url = \"", getSrcUrl(i),
                           "\" built = \"", getSrcBuilt(i), "\"/>", sep = "")
        writeLine(statement)
    }

#    load(file.path(.path.package("AnnBuilder"), "data", "AnnInfo.rda"))
    statement <- paste("<AnnBuilder:Entryid value = \"",
                       get(fileCol[1], AnnInfo)$short, "\"/>", sep = "")
    writeLine(statement)

    for(i in 2:length(fileCol)){
        statement <- paste("<AnnBuilder:Element value = \"",
                           fileCol[i], "\" describ = \"",
                           get(fileCol[i], AnnInfo)$short, "\"/>", sep = "")
        writeLine(statement)
    }

    writeLine("</AnnBuilder:Attr>")
}









