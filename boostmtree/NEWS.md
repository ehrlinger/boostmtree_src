# boostmtree 2.0.8

Cleveland Clinic Foundation patched build.

## New features

* `boostmtree()` and `partial.plot()` take `progress`, a function called as
  `progress(step, steps)` after each boosting iteration and after each partial-effect
  prediction. It reports progress where `verbose` cannot be seen: a knitted document
  captures the progress bar into its output, so a long render showed only which chunk
  it was in. An error in the callback becomes one warning and the fit carries on, and a callback cannot change the
  fit's random-number stream.

## Bug fixes

* `predict(..., partial = TRUE)` and `partial.plot()` no longer run out of
  memory on long time grids. `predict()` held every boosting step's linear
  predictor for every subject at every grid time at once, n.q x M x n vectors,
  in four such lists: the step path, its linear and mean transforms, and the
  extrapolated path in `boostmtree.predict.extrapolate()`. `partial.plot()`
  predicts every subject on every observed time, so an ordinal study with
  5,491 patients, 4,139 unique visit times and M = 100 asked for about 55 GB
  per list per grid value, and R died. Both paths now walk the boosting steps
  one at a time and keep a running sum, so a call holds about n.q x n vectors
  of the grid length. On a synthetic ordinal fit (250 subjects, 3,482 unique
  times, M = 30) peak memory for one partial prediction fell from 3.2 GB to
  0.4 GB. Results are unchanged, bit for bit. With test responses (`y`
  supplied), the path is walked twice, once for the error curve that picks
  the step and once for the fitted means at that step.
* `partial.plot()` predicts each grid value once and reads every requested
  response component from that prediction. It used to repeat the whole
  prediction for each component, three times over for a four-level ordinal
  response. It also skips the extrapolated `muhat` curves, which it never
  reads and which cost about as much again as the curves it does read.
  `predict()` itself still returns them.
* `predict()` matches each subject's times to the fitted time grid once per
  distinct time vector instead of once per subject. The nearest-time match is
  quadratic in the grid length, and under `partial = TRUE` every subject
  carries the whole grid, so it was about 90% of a partial prediction's time.
  With this and the changes above, the two-value synthetic ordinal
  `partial.plot()` fell from 87 s to 2.4 s, and one partial prediction from
  15.6 s to 2.2 s.
* Subjects on the same time grid share one time-design matrix in the returned
  predict object instead of each holding a copy.

# boostmtree 2.0.7

Cleveland Clinic Foundation patched build.

## Bug fixes

* Parallel cross-validated fits (`boostmtree.control(cv.cores > 1)`) no
  longer oversubscribe the host on a multithreaded BLAS such as OpenBLAS or
  MKL. 2.0.6 computes each terminal node's increment as one product over the
  node's stacked design rows, a product large enough for a threaded BLAS to
  start one thread per visible CPU in every forked worker: a 4-core fit at
  1000 subjects kept about 154 CPUs busy on a 192-CPU host and spent most of
  its time in the kernel. Each worker now runs its BLAS on one thread, and the
  parent drops to one thread while it forks and waits, then restores its
  previous count. Forking after a threaded BLAS has started its thread pool is
  also a known cause of spinning and hangs. A single-threaded BLAS, such as
  R's reference BLAS, is left alone.
* `RhpcBLASctl` joins `Imports`; it sets the BLAS thread count.

# boostmtree 2.0.6

Cleveland Clinic Foundation patched build.

## Performance

* Cross-validated fits on several cores use far less memory. Each subject's
  cross-validated view of every subject's linear predictor was a list of n
  short vectors, so the state held n^2 R objects per threshold. A garbage
  collection inside a forked worker touches every object, which turned the
  whole shared state into a private copy in each worker, and peak memory grew
  with the core count. The state is now one flat vector per subject, n objects
  in all, so a worker's collection no longer copies it, and each worker
  returns one vector per subject. The update adds one product per terminal
  node, over the node's stacked design rows, instead of one per subject. The
  arithmetic is unchanged: `mu`, `err.rate` and `m.opt` are identical to
  2.0.5. On the same number of cores a fit can take a few percent longer,
  because R collects garbage more often when it allocates long vectors.

## Tests

* The ordinal branch of the cross-validated update, which raises each
  threshold's state to the one below it, has a direct test, and an ordinal
  cross-validated fit is checked for identical results on 1 and 2 cores.

# boostmtree 2.0.5

Cleveland Clinic Foundation patched build.

## Performance

* Ordinal cross-validated fits spend far less time on one core. After the
  out-of-bag update, each threshold's state is raised to the one below it with
  `pmax()`. That ran once per (left-out, other) subject pair, n^2 calls per
  step, and `pmax()`'s fixed R-level cost per call made it the largest serial
  cost of the fit, about three quarters of the serial time at 600 subjects. It
  now runs once per subject over the subject's whole state. `pmax()` is
  elementwise, so the fit is unchanged: `mu`, `err.rate` and `m.opt` are
  identical to 2.0.4.

# boostmtree 2.0.4

Cleveland Clinic Foundation patched build.

## Performance

* `cv.flag = TRUE` fits can run the out-of-bag update on several cores with
  the new `boostmtree.control(cv.cores = )`. Each boosting step refits the
  terminal nodes once for every out-of-bag subject, which makes the
  cross-validated fit cost grow with the square of the number of subjects;
  those refits are independent within a step and are now forked with
  `parallel::mclapply()`. The default of 1 core keeps the sequential path.
  The refits draw no random numbers, so `mu`, `err.rate` and `m.opt` are
  identical for any core count. Windows, which cannot fork, uses 1 core.
* The out-of-bag update looks up each subject's terminal-node row once per
  step instead of once per (left-out, other) subject pair, and the node
  refit forms D'V^-1 once with `crossprod()` and reuses it for the Hessian
  and the score. Neither change alters the fit.
* `parallel` returns to `Imports`, now that it is used again.

# boostmtree 2.0.3

Cleveland Clinic Foundation patched build.

## Housekeeping

* Dropped the unused `parallel` dependency. The 2.0.0 refactor removed every
  `mclapply()` call, leaving `importFrom("parallel", "mclapply")` in
  `NAMESPACE` and `parallel` in `Imports` with nothing using them.
* CRAN Cookbook fixes: `boostmtree.news()` documents its return value as
  called for side effects and gains an example; the `Description` uses the
  lowercase `<doi:...>` form and fixes "is design for".

# boostmtree 2.0.2

Cleveland Clinic Foundation patched build.

## Bug fixes

* `na.action = "na.omit"` is now accepted. The formal is written
  `c("na.omit", "na.impute")[2]`, which evaluates to the single string
  `"na.impute"`, and a bare `match.arg()` takes its choices from that
  evaluated default. The documented value was therefore rejected, and the
  `na.omit` branch of the missing-data rule (drop subjects with any missing
  covariate, rather than only those missing all of them) was unreachable.
* `vimp.boostmtree(joint = TRUE)` no longer fails with "length of 'dimnames'
  [1] not equal to array extent". The collapse to a single `joint.vimp` label
  was applied only when `x.names` was supplied, so at the documented default
  the importance matrix had one row and one label per covariate. Supplying
  `x.names` took the working branch, which made the failure look intermittent.

Both are present in CRAN `boostmtree` 2.0.0 and have been reported upstream.

## Internal

* Removed `R/utilities.R`. The 2.0.0 refactor moved its contents into
  `boostmtree_math.R` and friends under `boostmtree.`-prefixed names and left
  the original file in the package; none of its 16 functions had a caller or
  was exported. It also carried an unguarded `solve()` that the live code path
  guards, which made it a standing source of false bug reports.
* Test coverage raised from 46.9% to 81.1%, with no source file left at zero.
  `partial.plot`, `marginal.plot`, `vimp.boostmtree`, `plot.boostmtree` and
  `print.boostmtree` had no tests at all; all plotting tests assert the absence
  of warnings. Both bug fixes above were found by writing these tests.

# boostmtree 2.0.1

Cleveland Clinic Foundation patched build.

## Bug fixes

* Fixed terminal-node coefficients (`gamma`) diverging when a model is fit
  with `cv.flag = TRUE`. The in-sample linear predictor and mean were
  refreshed only on the `cv.flag = FALSE` branch of the boosting loop, so
  with cross-validation enabled the working fit stayed frozen at its
  initialization. Residuals never shrank, the ridge penalty `lambda`
  collapsed, and each boosting step contributed a same-sized coefficient,
  so the linear predictor accumulated `M` unshrunk increments.
  Any output reading the stored `gamma` was affected: `predict()` with
  `use.cv.flag = FALSE`, `partial.plot()`, `marginal.plot()`,
  `predict(...)$muhat`, and `vimp()` on the non-CV path.
* Cross-validated output was not directly corrupted by the bug -- the
  cross-validation path accumulates its own out-of-bag increments -- but it
  was degraded indirectly, because the collapsed `lambda` also feeds
  `boostmtree.update.cv.step()`. Cross-validated results therefore move on
  this version, and the size and even the sign of the move depend on the
  configuration. Measured across 18 fits (`n` = 60 and 120; `M` = 50, 150
  and 250; three seeds each), the minimum of `err.rate` moved between 16.9%
  lower and 13.3% higher, and cross-validated `rmse` between 41.1% lower and
  15.6% higher; `lambda`, `phi` and `rho` shift as well. `m.opt` can move
  substantially -- it changed in 8 of the 18 fits, including 67 to 250 and
  119 to 248 -- so a stopping iteration chosen from an affected fit may have
  been badly truncated. Treat `err.rate`, `m.opt`, `rmse` and `mu` from any
  earlier `cv.flag = TRUE` fit as unusable: all of them must be re-derived
  from a re-fit on this version, and none may be carried over.
* Added a regression test asserting the in-sample boost converges under
  `cv.flag = TRUE`.

# boostmtree 2.0.0

This is a major maintenance and interface-cleanup release prepared for CRAN resubmission.
The release focuses on code correctness, documentation, naming consistency, and more predictable behavior for core fitting, prediction, and visualization workflows.

## Major changes

* Refactored the package around a modular source layout while preserving the core `ntree = 1` tree-based boosting path.
* Standardized the public interface to lower-case, dotted naming conventions.
  * Families are now specified as `"continuous"`, `"binary"`, `"nominal"`, and `"ordinal"`.
  * Returned object components now use the same naming convention.
* Added `boostmtree.control()` as the primary user-facing control interface for resampling, reproducibility, and other fitting options.
* Clarified first-pass scope: the current refactor supports `control$ntree = 1`.

## Fitting and model-object updates

* Reworked preprocessing so that subject-level covariates, identifiers, times, and responses remain aligned throughout fitting.
* Improved handling of non-continuous responses, including binary, nominal, and ordinal encoding.
* Improved support for univariate fits when `tm` and `id` are omitted.
* Restored the original tree-fitting logic more faithfully in key places while keeping the refactored package structure.
* Simplified and regularized the structure of returned objects.
  * Internal code keeps the boosted-subproblem (`q`) representation.
  * Public objects flatten single-response results consistently for continuous and binary fits.
* Removed legacy object components and hidden option paths that were no longer needed in the first-pass redesign.

## Resampling, OOB, and cross-validation

* Cleaned up the interaction between bootstrap sampling, OOB availability, cross-validation, and variable importance.
* `seed` now controls reproducibility only.
* Removed the silent legacy rewrite of `bootstrap = "none"` into a full in-bag user bootstrap.
* When `cv.flag = TRUE`, the fit now enforces an OOB-producing resampling rule and records OOB availability in the fitted object.
* Added explicit fitted-object fields documenting OOB behavior, including `oob.available` and `oob.subject.count`.

## Prediction, printing, and plotting

* Updated `predict.boostmtree()` to match the refactored object structure and naming conventions.
* Improved prediction handling for:
  * held-out longitudinal test data,
  * new subjects evaluated on the fitted training time grid,
  * user-supplied common time grids via partial prediction.
* Updated `print.boostmtree()` so grow and predict objects are summarized correctly and consistently.
* Updated `plot.boostmtree()` so plotting now goes to the active graphics device by default.
  * PDF output is still available when explicitly requested.
  * Plotting data can also be returned for user-directed graphics workflows.

## Variable importance and effect plots

* Redesigned variable-importance output around a classed `vimp.boostmtree` object.
* Replaced the old `vimpPlot()` workflow with `plot()` methods for variable-importance objects.
* Added clearer checks for grow-object variable importance when OOB information is unavailable.
* Reworked partial and marginal effect plotting around canonical function names:
  * `partial.plot()` / `partial.plot.boostmtree()`
  * `marginal.plot()` / `marginal.plot.boostmtree()`
* Updated effect-plot functions to use the active graphics device by default instead of creating PDF files automatically.
* Added response-label selection for partial and marginal plots in multi-level response settings.

## Documentation

* Rewrote the main `boostmtree` help page in a more user-facing style.
* Expanded examples to cover:
  * continuous longitudinal fits,
  * binary fits,
  * univariate fits,
  * held-out prediction,
  * AF-data illustration,
  * variable importance,
  * partial and marginal plots.
* Improved documentation for `print`, `plot`, `predict`, variable importance, and effect-plot methods.
* Clarified the interpretation of `phi`, `rho`, `lambda`, `mod.grad`, and prediction outputs such as `mu` and `muhat`.
* Updated package references to include later methodological and application papers:
  * Pande A., Ishwaran H., Blackstone E.H., Rajeswaran J., and Gillinov M. (2022). Application of gradient boosting in evaluating surgical ablation for atrial fibrillation. *SN Computer Science*, 3:466.
  * Pande A., Ishwaran H., and Blackstone E.H. (2022). Boosting for multivariate longitudinal responses. *SN Computer Science*, 3:186.

## Backward-incompatible changes

* The package now uses canonical lower-case, dotted names throughout the public interface.
* Returned object components were renamed to match the new naming convention.
* Legacy mixed-case and underscore-style argument names are no longer supported in the first-pass refactor.
* Legacy plotting helpers such as `vimpPlot()`, `partialPlot()`, and `marginalPlot()` were replaced by classed objects and `plot()` methods or lower-case dotted interfaces.
* Automatic PDF creation is no longer the default plotting behavior.

## Internal cleanup and bug fixes

* Removed conflicting and duplicated legacy helper code from `utilities.R`.
* Fixed multiple issues uncovered during example-driven testing, including bugs affecting:
  * univariate fits,
  * longitudinal prediction,
  * predict-object printing,
  * grow/predict consistency for `mu`, `muhat`, and related fields,
  * partial-plot indexing,
  * compatibility between stored `gamma` representations used during prediction.
* Removed legacy components no longer needed for the first-pass package design, including `forest.tol`.


# boostmtree 1.0.0

* Initial CRAN submission
