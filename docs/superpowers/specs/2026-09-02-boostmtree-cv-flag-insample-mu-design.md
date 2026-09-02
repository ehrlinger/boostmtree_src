# Fix: `cv.flag = TRUE` freezes the in-sample fit and corrupts `gamma`

Date: 2026-09-02
Repo: `ehrlinger/boostmtree_src` (fork of `kogalur/boostmtree`)
Branch: `fix/cv-flag-insample-mu`
Base: `8c5076c` (upstream v2.0.0, 2026-04-10)

## Problem

Fitting with `cv.flag = TRUE` produces a model whose stored terminal-node
coefficients (`gamma`) are unusable. Every downstream call that reads `gamma`
returns values far outside the observed response range:

- `predict.boostmtree()` with `use.cv.flag = FALSE` (the default for new data)
- `partial.plot()` and `marginal.plot()`
- `predict(...)$muhat`
- `vimp()` on the non-CV path

Observed severity: on simulated data with response range `[-28.9, 6.4]`,
`predict()` in-sample returns `[-167, 82]` at `M = 200`. The divergence is
linear in `M`, so it worsens the longer the model boosts.

### Root cause

`R/boostmtree_fit_tree.R:714` opens `if (model.info$cv.flag) {`, and `:741`
makes the in-sample branch its `else`. The refresh of the in-sample working
fit therefore runs *only* when `cv.flag = FALSE`:

```r
l.pred[[q]] <- lapply(..., function(i) l.pred.ref[[i]] + l.pred.db[[q]][[i]])
mu[[q]]     <- lapply(..., function(i) boostmtree.get.mu(l.pred[[q]][[i]], ...))
```

With `cv.flag = TRUE`, `mu[[q]]` stays frozen at its initialization for the
whole boosting loop. The consequences chain:

1. `boostmtree.estimate.lambda()` reads the frozen `mu` (via `mu.for.lambda`,
   which is the in-sample `mu` under the default `control$cv.lambda = FALSE`)
   and so sees full-size residuals at every iteration.
2. `boostmtree.update.gls.parameters()` likewise inflates `phi`.
3. The ridge penalty `lambda` collapses by ~3 orders of magnitude.
4. With effectively no penalty, each terminal-node coefficient vector is fit
   to the full response rather than to a shrinking residual, so `max|gamma|`
   never decays.
5. `l.pred.db` — which *is* updated each iteration — accumulates `M` unshrunk
   increments, producing a random walk instead of a converging boost.

Measured, same data and seed, `M = 200`:

| | `cv.flag=FALSE` | `cv.flag=TRUE` |
|---|---|---|
| `phi` | 0.667 | 2.303 |
| `lambda` | ~60,000 | 39 |
| `max\|gamma\|` at m = 1, 20, 100, 200 | 3.75, 1.48, 0.081, 0.040 | 3.75, 3.30, 3.41, 3.35 |
| `predict()` in-sample range | inside observed y | `[-167, 82]` |

### Why it went unnoticed

When `cv.flag = TRUE`, `boostmtree.finalize.fit()` overwrites `fit.info$mu`
at `:410` with the *cross-validated* `mu`. The CV path accumulates only each
subject's out-of-bag steps via `gamma.i.list`, and empirically survives the
collapsed `lambda`: CV error curves differ by <0.1% between the broken and
fixed builds. So `fit$mu`, `fit$err.rate`, `fit$m.opt` and CV RMSE all look
correct — the bug is invisible until `gamma` is read.

`partial.plot()` compounds this: `R/partial_plot.boostmtree.R:359` silently
downgrades `use.cv.flag` to `FALSE` whenever `subset` is non-`NULL`, so
passing a subset forces callers onto the corrupt path. Unsubsetted partial
plots average the divergence toward the mean and look plausible; subsetted
ones do not.

## Scope

`R/boostmtree_fit_tree.R:714` is the only defective gate. Audited siblings:

- `:379` — `finalize.fit` returning in-sample `mu` when `cv.flag = FALSE`. Correct.
- `:518` — `cv.state` initialization. Correct.
- `:694` — CV-only `gamma.i.list` update. Correct.
- `:781` — OOB-count warning. Correct.

## Design

### Change

Make the in-sample update unconditional. The CV state is *additional*
bookkeeping; the in-sample linear predictor must advance every iteration
because `lambda`, `rho` and `phi` estimation all read it.

```r
if (model.info$cv.flag) {
  ...CV state update...          # unchanged
}
...in-sample l.pred / mu update...   # was the `else` body; now always runs
```

Mechanically: `:741` `} else {` becomes `}`, and the matching brace at `:757`
is removed. No logic is rewritten and no control flow is added.

Cost under `cv.flag = TRUE`: one additional `lapply` over subjects per
boosting iteration, negligible against the per-iteration `rfsrc` fit.

Behaviour deliberately left unchanged: `finalize.fit` still overwrites
`fit$mu` with the CV `mu` when `cv.flag = TRUE`. That is intended.

### Version

`2.0.0` -> `2.0.1` in `boostmtree/DESCRIPTION`, with a matching `NEWS.md`
entry. Patch digit only, straight three-digit. This gives the analysis and
its methods statement an identifier distinct from upstream 2.0.0.

A non-numeric version (`2.0.0.fix`, `2.0.0.ccf`) is not an option: R's
`package_version()` rejects letters outright and the build would fail to
install. A fourth numeric component (`2.0.0.9000`) parses but violates the
project's straight-three-digit rule.

Fork identity is instead carried where it is legal and durable:

- git tag `v2.0.1-ccf` (tags permit letters), so the installable ref is marked
- a sentence in the DESCRIPTION `Description:` prose noting the CCF patch
- the `Remote*` fields that `remotes` stamps into the installed DESCRIPTION

That last one is decisive for provenance: after install,
`packageDescription("boostmtree")$RemoteSha` pins the exact commit and
`$RemoteUsername` returns `ehrlinger`, while an upstream CRAN 2.0.1 would
carry no `Remote*` fields at all. Cite the SHA in the methods statement.

```r
remotes::install_github("ehrlinger/boostmtree_src",
                        subdir = "boostmtree", ref = "v2.0.1-ccf")
```

The package lives in the `boostmtree/` subdirectory and the repo has no root
`DESCRIPTION`, so repo-root files (including this spec) are outside the build
and cannot affect `R CMD check` or the tarball. No restructuring required.

### Testing

The repo has no `tests/` directory; this adds `testthat` scaffolding:

- `boostmtree/tests/testthat.R`
- `boostmtree/tests/testthat/test-cv-flag-gamma-decay.R`
- `Suggests: testthat (>= 3.0.0)` and `Config/testthat/edition: 3` in `DESCRIPTION`

The regression test is the one that would have caught this: fit a small
model with `cv.flag = TRUE` and assert that the boosting actually converges —

1. `max|gamma|` at the final `m` is < 40% of its value at `m = 1`.
   Measured at the test's own settings (n=60, N=5, M=100): 0.204 when fixed
   versus 1.133 when broken. The 0.40 threshold sits with >2x margin on both
   sides, so it separates cleanly without being brittle to RNG or platform.
2. `predict(fit, use.cv.flag = FALSE)$mu` lies within the observed response
   range expanded by 25% of its width on each side.

Asserting on `gamma` decay rather than on a stored numeric fixture keeps the
test robust to RNG and platform differences while still failing loudly if the
in-sample path is ever frozen again.

Kept small (low `n`, `N`, `M`) to stay inside a sane check-time budget.

## Verification gate

All must pass before the branch is considered done:

1. **No collateral change.** `cv.flag = FALSE` results bit-identical to
   upstream `8c5076c`. This is the path any fallback analysis would use.
2. **Fix works.** Under `cv.flag = TRUE`: `max|gamma|` decays across `m`,
   `lambda` recovers to its penalized scale, and `predict()` in-sample lands
   inside the observed response range.
3. **Package health.** `R CMD check` shows no new WARNINGs or NOTEs against
   an upstream baseline run.

## Out of scope

- The hvtiPlotR subgroup ensemble-curve helper (follow-on; both outcomes are
  separate univariate fits, so `n.q = 1` and only the continuous path applies).
- The `partial.plot()` `subset` -> `use.cv.flag = FALSE` downgrade at
  `partial_plot.boostmtree.R:359`. It is defensible on its own terms — CV
  coefficients are subject-indexed and genuinely do not survive subsetting.
  Worth reporting upstream, but not changed here.
- Upstream issue and PR text against `kogalur/boostmtree`.
