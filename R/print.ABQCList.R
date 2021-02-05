# This function prints the quality control results for a givan data package.
#
# Copyright 2002, J. Zhang. All rights reserved.
#

print.ABQCList <- function(x, ...){
    cat("\n",
        paste("\nQuality control information for ",
              x$name, ":\n", sep = ""),
        paste("\nDate built:", x$built),
        ifelse(is.null(x$probeNum), "",
        paste("\nNumber of probes:", x$probeNum)),
        ifelse(is.null(x$numMissMatch), "",
        paste("\nProbe number missmatch:", paste(x$numMissMatch,
              sep = "", collapse = "; "))),
        ifelse(is.null(x$probeMissMatch), "",
        paste("\nProbe missmatch:", paste(x$probeMissMatch,
              sep = "", collapse = "; "))),
        ifelse(is.null(x$probeMapped), "", paste(
        "\nMappings found for probe based rda files: \n",
        paste("\t", names(x$probeMapped), "found", x$probeMapped,
                    "of", x$probeNum, sep = " ", collapse = "\n"),
                                   sep = "", collapse ="\n")),
        "\nMappings found for non-probe based rda files:\n",
        paste("\t", names(x$otherMapped), "found", x$otherMapped,
              sep = " ", collapse = "\n"),"\n\n")
}

