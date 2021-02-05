pfamAC2SCOP <- function(ac=NULL){
    getPfam(pfamCon=getPfamDb(), FROM="AC", TO=c("SCOP", "PLACEMENT"), TABLE="SCOP", FOR=ac)
}
pfamSCOP2AC <- function(scop=NULL){
    getPfam(pfamCon=getPfamDb(), FROM="SCOP", TO="AC", TABLE="SCOP", FOR=scop)
}

pfamAC2PDB <- function(ac=NULL){
    getPfam(pfamCon=getPfamDb(), FROM="AC", TO=c("PDB", "startPoint", "endPoint"), TABLE="PDB", FOR=ac)
}
pfamPDB2AC <- function(pdb=NULL){
    getPfam(pfamCon=getPfamDb(), FROM="PDB", TO="AC", TABLE="PDB", FOR=pdb)
}
