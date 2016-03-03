## #########################################################################
##
##  Generic Header
##
## #########################################################################

options(object.size=Inf,expressions=100000,memory=Inf,width=150)

output.path = paste(getwd(), "/output/", sep="")

if (!file.exists(path = output.path)) {
  dir.create(path = output.path, showWarnings = TRUE)
}

use.package = gsub("-use.package=","",grep("-use.package=",commandArgs(),value=T));

if ( length(use.package) == 0) {
  use.package = FALSE
}
if (use.package == "TRUE") {
  use.package = TRUE
} else {
  use.package = FALSE
}

if (use.package) {
  library("_PROJECT_PACKAGE_NAME_")
} else {
  source("load.R")
}

update.benchmark = gsub("-update.benchmark=","",grep("-update.benchmark=",commandArgs(),value=T));

if ( length(update.benchmark) == 0) {
  update.benchmark = FALSE
}
if (update.benchmark == "TRUE") {
  update.benchmark = TRUE
} else {
  update.benchmark = FALSE
}

stealth.flag = gsub("-stealth.flag=","",grep("-stealth.flag=",commandArgs(),value=T));

if ( length(stealth.flag) == 0) {
  stealth.flag = FALSE
}
if (stealth.flag == "TRUE") {
  stealth.flag = TRUE
} else {
  stealth.flag = FALSE
}

cat("\n use.package = ",      use.package, "\n")
cat("\n update.benchmark = ", update.benchmark, "\n")
cat("\n stealth.flag = ",     stealth.flag, "\n")

