pubDataURLs <- read.table("pubDataURLs.txt", colClasses=rep("character", 2),
                          header=TRUE)
pubDataURLs <- lapply(list(pubDataURLs), function(x) {
    urls = x$url
    names(urls) = x$name
    as.list(urls)
})
pubDataURLs <- pubDataURLs[[1]]

