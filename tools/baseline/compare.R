# Fails loudly if the cv.flag = FALSE path moved at all.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript compare.R <baseline.rds> <patched.rds>")

baseline <- readRDS(args[1])
patched <- readRDS(args[2])

components <- names(baseline)
mismatched <- character(0)

for (nm in components) {
  if (!identical(baseline[[nm]], patched[[nm]])) {
    mismatched <- c(mismatched, nm)
  }
}

if (length(mismatched) > 0L) {
  cat("MISMATCH in:", paste(mismatched, collapse = ", "), "\n")
  quit(status = 1L)
}

cat("cv.flag = FALSE path bit-identical across", length(components),
    "components:", paste(components, collapse = ", "), "\n")
