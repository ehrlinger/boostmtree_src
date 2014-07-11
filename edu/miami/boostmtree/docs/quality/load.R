###########################################################################
#
#  Load Working Copies of files into an R session.
###########################################################################

###########################################################################
#  Check for valid working directory. This file must
#  be sourced from the build/regression directory.
###########################################################################
build.path = NULL

if (file.exists(path = paste(getwd(), "/../quality", sep="")) &
    file.exists(path = paste(getwd(), "/R", sep="")) &
    file.exists(path = paste(getwd(), "/lib", sep="")) &
    file.exists(path = paste(getwd(), "/../../build", sep=""))) {

  build.path = paste(getwd(), "/../../", sep="")
  argh.path = paste(getwd(), "/R/", sep="")
  lib.path = paste(getwd(), "/lib/_PROJECT_ARCH_TYPE_/", sep="")
}

if (is.null(build.path)) {
  stop("\nBMT ERROR:  Problem loading benchmark files.  Working Directory Invalid:  ",
       "\n", getwd(),
       "\nPlease run from build/quality.  ")
}

## Dependency automatically loaded when package is loaded.
library("randomForestSRCM")

###########################################################################
#  Source files, and load native code object.
###########################################################################

source(paste(argh.path, "boostmtree.R", sep=""))
source(paste(argh.path, "predict.boostmtree.R", sep=""))
source(paste(argh.path, "generic.predict.boostmtree.R", sep=""))
source(paste(argh.path, "partialPlot.R", sep=""))
source(paste(argh.path, "plot.boostmtree.R", sep=""))
source(paste(argh.path, "print.boostmtree.R", sep=""))
source(paste(argh.path, "utilities.R", sep=""))


source(paste(argh.path, "boostmtree.news.R", sep=""))
source(paste(argh.path, "zzz.R", sep=""))

## dyn.load(paste(lib.path, "_PROJECT_LIBRARY_NAME_", sep=""))

cat("\n           Version:  _PROJECT_VERSION_ID_ ");
cat("\n Development Build:  _PROJECT_BUILD_ID_ ");
cat("\n      Architecture: _PROJECT_ARCH_TYPE_ ");
cat("\n   Affirm Parallel: _PROJECT_SUPPORT_OPENMP_ \n\n");
