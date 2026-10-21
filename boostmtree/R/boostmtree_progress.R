# A caller's progress callback, called as progress(step, steps). A callback that
# fails warns once and is not called again, so reporting can never stop a fit.
boostmtree.progress.caller <- function(progress) {
  if (is.null(progress)) return(NULL)
  if (!is.function(progress)) stop("`progress` must be NULL or a function(step, steps).")
  state <- new.env(parent = emptyenv())
  state$live <- TRUE
  function(step, steps) {
    if (!state$live) return(invisible(NULL))
    # a callback must not change the fit, so whatever random numbers it draws are put back
    seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", seed, envir = globalenv()))
    tryCatch(progress(step, steps), error = function(e) {
      state$live <- FALSE
      warning("`progress` failed and will not be called again: ", conditionMessage(e), call. = FALSE)
    })
    invisible(NULL)
  }
}
