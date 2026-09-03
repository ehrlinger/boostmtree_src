# boostmtree: Cleveland Clinic Foundation patched fork

## What this is

This is a fork of [`kogalur/boostmtree`](https://github.com/kogalur/boostmtree)
carrying one correctness fix. In upstream `boostmtree`, fitting with
`cv.flag = TRUE` skipped the in-sample `l.pred`/`mu` refresh inside the boosting
loop. The in-sample residuals therefore never shrank: `phi` inflated, the ridge
penalty `lambda` collapsed by roughly three orders of magnitude, and every
boosting step contributed a same-sized terminal-node coefficient, so the linear
predictor accumulated `M` unshrunk increments and diverged linearly in `M`.
Anything that reads the non-cross-validated path (`predict(..., use.cv.flag =
FALSE)`, `partial.plot()`, `vimp()`) was affected. The `cv.flag = FALSE` path
was never wrong and is bit-identical before and after the fix. The package
itself lives in the [`boostmtree/`](boostmtree) subdirectory of this repository.

## Install

```r
remotes::install_github("ehrlinger/boostmtree_src",
                        subdir = "boostmtree",
                        ref    = "v2.0.1-ccf")
```

If the `v2.0.1-ccf` tag is not yet published, install from the branch instead,
or from a local clone with `remotes::install_local("boostmtree")` run from the
repository root. Check the
[tags page](https://github.com/ehrlinger/boostmtree_src/tags) for what is
currently available.

## Verify your machine reproduces the reference fit

This matters. `boostmtree` obtains terminal-node memberships from
`randomForestSRC::predict.rfsrc(..., ptn.count = k)$ptn.membership`.
`DESCRIPTION` only requires `randomForestSRC (>= 3.5.0)`, but every verification
run for this fork was done against **randomForestSRC 3.6.2, R 4.6.1, macOS
arm64**. A different `randomForestSRC` can partition the covariate space
differently, which changes `gamma`, which changes the fitted values, silently,
with no error raised. Run the check before you trust any numbers:

```sh
Rscript tools/baseline/reproducibility-check.R
```

It refits the committed fixed-seed reference model, compares all five
fingerprint components (`mu`, `lambda`, `gamma`, `phi`, `rho`) against
`tools/baseline/baseline-cvflag-false.rds` at `1e-8` relative tolerance, prints
your R version, platform, `boostmtree`/`randomForestSRC`/`nlme` versions and
BLAS, and exits non-zero on any divergence. The tolerance is deliberately loose
enough to absorb legitimate last-bit differences from a different BLAS/LAPACK or
CPU, and far tighter than the O(1) shift a different tree partition would cause.

**If it fails:** your environment does not reproduce the reference fit, and
fitted values from this installation may differ from the verified results. Do
not proceed on those numbers. Compare the printed environment against the
reference, pin `randomForestSRC` to the reference version
(`remotes::install_version("randomForestSRC", version = "3.6.2")`), and re-run.
If it still fails with `randomForestSRC` pinned, report the full output to the
package maintainer rather than using the build.

(`tools/baseline/compare.R` is a different tool: it uses `identical()` and is for
proving a code change moved nothing *on one machine*. It is not appropriate
across machines.)

## Performance: use `cv.flag = TRUE` once, then turn it off

`cv.flag = TRUE` is about **5x slower** than `cv.flag = FALSE`, measured at
10.1 s versus 1.9 s for `M = 50`, `n = 150`. Recommended workflow:

1. Fit once with `cv.flag = TRUE` to select the optimal number of boosting
   steps `M` (`fit$m.opt`).
2. Refit with `cv.flag = FALSE` at that `M`, and use *that* fit for all
   downstream `predict()`, `partial.plot()` and `vimp()` work.

```r
cvfit <- boostmtree(x, tm, id, y, family = "continuous", M = 500, cv.flag = TRUE)
fit   <- boostmtree(x, tm, id, y, family = "continuous", M = cvfit$m.opt,
                    cv.flag = FALSE)
```

## Warning: re-derive results from any pre-fix `cv.flag = TRUE` fit

`err.rate`, `m.opt`, `rmse` and `mu` produced by a `cv.flag = TRUE` fit from an
**unpatched** `boostmtree` are affected by the bug and must be **re-derived with
this build, not carried over**. That includes any saved `.rds` fit objects,
cached model selections, and any `M` chosen from a pre-fix cross-validation
curve. Refit from the data.

## Continuous integration

* `.github/workflows/R-CMD-check.yaml`: `R CMD check --as-cran` on
  ubuntu/macos/windows, R release, including the `testthat` suite.
* `.github/workflows/install-smoke.yaml`: installs the package from the
  subdirectory the way a consumer would and runs `tools/smoke/install-smoke.R`,
  a small fixed-seed `cv.flag = TRUE` fit whose `predict(use.cv.flag = FALSE)`
  output must stay inside a sane band around the observed response range.

Both live at the repository root and are outside the package build, so nothing
here ships in the tarball.
