## #########################################################################
##
##  Script to compress the data files distributed in the package.
##  Used during CRAN package preparation.
##
## #########################################################################

library(tools)

fileSuffix <- c("txt", "tab", "csv")

for (j in 1:length(fileSuffix)) {

    fileName <- dir(path=".", pattern=fileSuffix[j])
    filePrefix <- sapply(fileName, function(x) strsplit(x, paste(".", fileSuffix[j], sep="")))

    dataList <- sapply(fileName, function(x) read.table(x, header = TRUE))

    if (length(dataList) > 0) {
        for (i in 1:length(dataList)) {
            assign(as.character(filePrefix[i]), dataList[[i]])
            save(list=as.character(filePrefix[i]), file=paste(filePrefix[i], ".rda", sep=""))
        }
    }
}

resaveRdaFiles(".")
