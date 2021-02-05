rejectORF <- function(pkgPath, pkgName="YEAST", orfDecisionUrl="http://www.broad.mit.edu/ftp/pub/annotation/fungi/comp_yeasts/S6.RFC_test/a.orf_decisions.txt"){
    orfUrl <- url(orfDecisionUrl)
    orf <- read.delim(orfUrl, as.is=TRUE)
    YEASTREJECTORF <- orf$name[orf$RFC=="reject"]
    save(YEASTREJECTORF, file=file.path(pkgPath, pkgName, "data", "YEASTREJECTORF.rda"), compress=TRUE)
}
