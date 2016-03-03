source("header.R")


lsp <- function(package, all.names = FALSE, pattern) 
{
  package <- deparse(substitute(package))
  ls(
      pos = paste("package", package, sep = ":"), 
      all.names = all.names, 
      pattern = pattern
  )
}


run.example <- function(topic) {
  eval(parse(text=paste("example(", topic, ", echo = FALSE, ask = FALSE, run.dontrun = TRUE)", sep="")))
}

sink(paste(output.path, "all.examples.txt", sep=""), append=FALSE)
sink()

if (TRUE) {

  all.names <- lsp(_PROJECT_PACKAGE_NAME_)

  for (i in 1:length(lsp)) {
    
    sink(paste(output.path, "all.examples.txt", sep=""), append=TRUE)
    cat("Executing", all.names[i], "() ...", "\n\n")
    sink()
    run.example(eval(all.names[i]))
    cat("\n\n\n\n\n")
  
  }
}

sink(paste(output.path, "all.examples.txt", sep=""), append=TRUE)
print(TRUE)
cat ("  :  returned nominally \n\n")
cat ("\n\n")
sink()

