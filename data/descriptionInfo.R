descriptionInfo <- read.delim("descriptionInfo.txt", as.is=TRUE, header=TRUE, row.names=NULL)

biocViewsIndex <- which(colnames(descriptionInfo)=="biocViews")

organismIndex <- which(colnames(descriptionInfo)=="organism")
descriptionInfo[, biocViewsIndex] <-   apply(descriptionInfo, 1,
                                               function(pkgInfo){
                                                   biocViewsOrgName <- gsub(" ", "_", pkgInfo[organismIndex])
                                                   if(biocViewsOrgName != ""){
                                                       paste(pkgInfo[biocViewsIndex], biocViewsOrgName, sep=", ")
                                                   }else{
                                                       return(pkgInfo[biocViewsIndex])
                                                   }
                                               })


manufacturerIndex <- which(colnames(descriptionInfo)=="manufacturer")
descriptionInfo[, biocViewsIndex] <- apply(descriptionInfo, 1,
                                           function(pkgInfo){
                                               if(pkgInfo[manufacturerIndex] != ""){
                                                   biocViewsManuName <- paste(pkgInfo[manufacturerIndex], "Chip", sep="")
                                                   paste(pkgInfo[biocViewsIndex], biocViewsManuName, sep=", ")
                                               }else{
                                                   return(pkgInfo[biocViewsIndex])
                                               }
                                           })


biocPkgNameIndex <- which(colnames(descriptionInfo)=="biocPkgName")
descriptionInfo[, biocViewsIndex] <- apply(descriptionInfo, 1,
                                           function(pkgInfo){
                                               notAnnPkgKey <- c("probe$", "cdf$", "CHRLOC$", "homology$", "^cMAP$", "^GO$", "^KEGG$", "^PFAM$", "^YEAST$", "^$")
                                               numMatchKey <- sum(unlist(lapply(notAnnPkgKey, function(key) grep(key, pkgInfo[biocPkgNameIndex]))))
                                               if(numMatchKey==0){
                                                   paste(pkgInfo[biocViewsIndex], pkgInfo[biocPkgNameIndex], sep=", ")
                                               }else{
                                                   return(pkgInfo[biocViewsIndex])
                                               }
                                           })


