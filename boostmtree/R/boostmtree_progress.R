# A caller's progress callback, called as progress(step, steps). A callback that
# fails warns once and is not called again, so reporting can never stop a fit.
boostmtree.progress.caller <- function(progress) {
  if (is.null(progress)) return(NULL)
  if (!is.function(progress)) stop("`progress` must be NULL or a function(step, steps).")
  state <- new.env(parent = emptyenv())
  state$live <- TRUE
  function(step, steps) {
    if (!state$live) return(invisible(NULL))
    tryCatch(progress(step, steps), error = function(e) {
      state$live <- FALSE
      warning("`progress` failed and will not be called again: ", conditionMessage(e), call. = FALSE)
    })
    invisible(NULL)
  }
}
