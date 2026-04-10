predict.boostmtree <- function(
  object,
  x,
  tm,
  id,
  y,
  M = NULL,
  eps = 1e-5,
  use.cv.flag = FALSE,
  partial = FALSE,
  ...
) {
  call <- match.call(expand.dots = TRUE)
  helper <- if (exists("generic.predict.boostmtree", mode = "function", inherits = TRUE)) {
    get("generic.predict.boostmtree", mode = "function", inherits = TRUE)
  } else {
    ns <- tryCatch(asNamespace("boostmtree"), error = function(e) NULL)
    if (!is.null(ns) && exists("generic.predict.boostmtree", envir = ns, inherits = FALSE)) {
      get("generic.predict.boostmtree", envir = ns, inherits = FALSE)
    } else {
      stop("Internal prediction helper `generic.predict.boostmtree()` is not available.")
    }
  }
  call[[1L]] <- helper
  eval(call, envir = parent.frame())
}
