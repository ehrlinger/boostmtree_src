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
