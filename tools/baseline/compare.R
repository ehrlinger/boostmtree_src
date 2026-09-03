# Fails loudly if the cv.flag = FALSE path moved at all.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript compare.R <baseline.rds> <patched.rds>")

baseline <- readRDS(args[1])
patched <- readRDS(args[2])

components <- names(baseline)

# Guard against a vacuous pass: an empty fingerprint, or one whose components
# are all NULL (e.g. regenerated from a univariate fit, where lambda/phi/rho
# are absent), would otherwise compare identical(NULL, NULL) and report success.
if (length(components) == 0L) {
  cat("ERROR: baseline fingerprint has no components; nothing was compared.\n")
  quit(status = 1L)
}
null.components <- components[vapply(components,
                                     function(nm) is.null(baseline[[nm]]),
                                     logical(1))]
if (length(null.components) > 0L) {
  cat("ERROR: baseline component(s) NULL, comparison would be vacuous:",
      paste(null.components, collapse = ", "), "\n")
  quit(status = 1L)
}

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
