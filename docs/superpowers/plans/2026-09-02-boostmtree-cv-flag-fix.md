# boostmtree `cv.flag` In-Sample Mu Fix — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop `cv.flag = TRUE` from freezing the in-sample fit, so the stored `gamma` — and therefore `predict()`, `partial.plot()`, `marginal.plot()`, `muhat` and `vimp()` — returns values inside the observed response range.

**Architecture:** A three-line control-flow correction in the boosting loop of `R/boostmtree_fit_tree.R`. The `else` that gates the in-sample `l.pred`/`mu` refresh is removed so the refresh always runs; the CV-state block above it is untouched. Guarded by a new testthat regression test that asserts the boosting actually converges, and by a byte-level fingerprint proving the `cv.flag = FALSE` path is unchanged.

**Tech Stack:** R (>= 4.3.0), randomForestSRC (>= 3.5.0), testthat 3e, `remotes` for install.

## Global Constraints

- Repo: `/Users/ehrlinj/Documents/GitHub/boostmtree_src`, branch `fix/cv-flag-insample-mu`, base `8c5076c`.
- The R package lives in the `boostmtree/` **subdirectory**. Repo-root paths (`docs/`, `tools/`) are outside the package build — never move package files to the root.
- Version becomes exactly `2.0.1`. Straight three digits. No fourth component, no letters — `package_version()` rejects letters and the build would fail to install.
- Fork identity is carried by the git tag `v2.0.1-ccf`, a sentence in the DESCRIPTION `Description:` prose, and the `Remote*` fields `remotes` stamps at install. Not by the version string.
- Do not modify `R/partial_plot.boostmtree.R`. The `subset` -> `use.cv.flag = FALSE` downgrade at `:359` is out of scope.
- Do not change `boostmtree.finalize.fit()` at `:379`/`:410`. Overwriting `fit$mu` with the CV mu under `cv.flag = TRUE` is intended behaviour.
- Every commit message ends with the trailer `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
- Never push to `master` on `kogalur/boostmtree`. All work stays in the `ehrlinger/boostmtree_src` fork.

## Scratch library convention

Tasks install the package into a throwaway library so the user's main R
library is never mutated. Export this in every shell that runs R:

```bash
export BMLIB=/private/tmp/claude-504/-Users-ehrlinj-Documents-GitHub-boostmtree/5afbe653-7bf3-476e-b318-a4c4d253359b/scratchpad/bmlib
mkdir -p "$BMLIB"
```

Install with:

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R CMD INSTALL --no-docs --no-byte-compile -l "$BMLIB" boostmtree
```

Run R with that library first on the path via `R_LIBS_USER="$BMLIB"`.

---

### Task 1: Capture the upstream baseline

Must run **before** any source change — it fingerprints unpatched behaviour so Task 4 can prove the `cv.flag = FALSE` path did not move.

**Files:**
- Create: `tools/baseline/fingerprint.R`
- Create: `tools/baseline/baseline-cvflag-false.rds` (generated)
- Create: `tools/baseline/upstream-check.log` (generated)

**Interfaces:**
- Consumes: nothing.
- Produces: `tools/baseline/fingerprint.R`, runnable as
  `Rscript tools/baseline/fingerprint.R <output.rds>`; writes a list with
  elements `mu`, `lambda`, `gamma`, `phi`, `rho` from a fixed-seed
  `cv.flag = FALSE` fit. Task 4 reads this file.

- [ ] **Step 1: Write the fingerprint script**

Create `tools/baseline/fingerprint.R`:

```r
# Deterministic fingerprint of the cv.flag = FALSE fitting path.
# Used to prove the cv.flag fix introduces no collateral change.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript fingerprint.R <output.rds>")

suppressMessages(library(boostmtree))

set.seed(20260902)
d <- simLong(
  n = 60, n.time = 5, rho = 0.8, model = 2,
  family = "continuous", q = 0
)$data.list

fit <- boostmtree(
  d$features, d$time, d$id, d$y,
  family = "continuous", M = 50, cv.flag = FALSE
)

saveRDS(
  list(mu = fit$mu, lambda = fit$lambda, gamma = fit$gamma,
       phi = fit$phi, rho = fit$rho),
  args[1]
)
cat("fingerprint written to", args[1], "\n")
```

- [ ] **Step 2: Install upstream code into the scratch library**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
export BMLIB=/private/tmp/claude-504/-Users-ehrlinj-Documents-GitHub-boostmtree/5afbe653-7bf3-476e-b318-a4c4d253359b/scratchpad/bmlib
mkdir -p "$BMLIB"
R CMD INSTALL --no-docs --no-byte-compile -l "$BMLIB" boostmtree
```

Expected: ends with `* DONE (boostmtree)`.

- [ ] **Step 3: Generate the baseline fingerprint**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R_LIBS_USER="$BMLIB" Rscript tools/baseline/fingerprint.R tools/baseline/baseline-cvflag-false.rds
```

Expected: `fingerprint written to tools/baseline/baseline-cvflag-false.rds`

- [ ] **Step 4: Capture the upstream `R CMD check` baseline**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R_LIBS_USER="$BMLIB" R CMD check --no-manual --no-build-vignettes boostmtree \
  > tools/baseline/upstream-check.log 2>&1
grep -E "^(\*|Status)" tools/baseline/upstream-check.log | tail -20
rm -rf boostmtree.Rcheck
```

Record the final `Status:` line. Task 6 compares against it. A non-clean
upstream baseline is fine and expected — the gate is "no *new* problems",
not "zero problems".

- [ ] **Step 5: Commit**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
git add tools/baseline
git commit -m "test: capture upstream baseline before cv.flag fix

Fingerprints the cv.flag = FALSE path and records the upstream
R CMD check status, so the fix can be proven free of collateral change.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Add testthat scaffolding and the failing regression test

**Files:**
- Modify: `boostmtree/DESCRIPTION` (add `Suggests`, `Config/testthat/edition` after line 11)
- Create: `boostmtree/tests/testthat.R`
- Create: `boostmtree/tests/testthat/test-cv-flag-gamma-decay.R`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: a test file that MUST fail against unpatched source and pass after Task 3.

- [ ] **Step 1: Add testthat to DESCRIPTION**

Insert after line 11 (`Imports: ...`), before line 12 (`Description: ...`):

```
Suggests: testthat (>= 3.0.0)
Config/testthat/edition: 3
```

- [ ] **Step 2: Create the testthat runner**

Create `boostmtree/tests/testthat.R`:

```r
library(testthat)
library(boostmtree)

test_check("boostmtree")
```

- [ ] **Step 3: Write the failing regression test**

Create `boostmtree/tests/testthat/test-cv-flag-gamma-decay.R`:

```r
# Regression test for the cv.flag = TRUE in-sample mu freeze.
#
# When the in-sample l.pred/mu refresh is skipped under cv.flag = TRUE, the
# residuals never shrink, phi inflates, the ridge penalty lambda collapses by
# ~3 orders of magnitude, and every boosting step contributes a same-sized
# gamma. The linear predictor then accumulates M unshrunk increments and
# diverges linearly in M.
#
# Both assertions below test that the boosting actually converges, rather than
# pinning stored numerics, so they stay robust across RNG and platform.

test_that("cv.flag = TRUE still yields a converging in-sample boost", {
  set.seed(202609)
  d <- simLong(
    n = 60, n.time = 5, rho = 0.8, model = 2,
    family = "continuous", q = 0
  )$data.list

  fit <- boostmtree(
    d$features, d$time, d$id, d$y,
    family = "continuous", M = 100, cv.flag = TRUE
  )

  # 1. Terminal-node coefficients must decay as the residual shrinks.
  #    Measured: 0.204 when correct, 1.133 when frozen.
  gamma.q <- fit$gamma[[1]]
  node.max <- function(m) max(abs(as.matrix(gamma.q[[m]])[, -1, drop = FALSE]))
  decay.ratio <- node.max(length(gamma.q)) / node.max(1)

  expect_lt(decay.ratio, 0.40)

  # 2. In-sample predictions must stay near the observed response range.
  #    The non-CV path is what predict(), partial.plot() and vimp() read.
  observed <- range(d$y)
  width <- diff(observed)
  lower <- observed[1] - 0.25 * width
  upper <- observed[2] + 0.25 * width

  predicted.mu <- unlist(predict(fit, use.cv.flag = FALSE)$mu)
  expect_true(length(predicted.mu) > 0)
  expect_true(all(is.finite(predicted.mu)))

  predicted <- range(predicted.mu)

  expect_gte(predicted[1], lower)
  expect_lte(predicted[2], upper)
})
```

- [ ] **Step 4: Run the test and verify it FAILS**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R CMD INSTALL --no-docs --no-byte-compile -l "$BMLIB" boostmtree
R_LIBS_USER="$BMLIB" Rscript -e 'testthat::test_local("boostmtree", reporter="summary")'
```

Expected: FAIL. `decay.ratio` is ~1.13, not < 0.40, and the predicted range
(~`[-58.8, 38.9]`) falls outside the allowed band (~`[-31.2, 13.4]`).
Both `expect_lt` and the range expectations should report failures.

If the test PASSES here, stop — the working tree is not at the unpatched
baseline. Check `git log --oneline -1` and `git status`.

- [ ] **Step 5: Commit the failing test**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
git add boostmtree/DESCRIPTION boostmtree/tests
git commit -m "test: add failing regression test for cv.flag gamma decay

Asserts the in-sample boost converges under cv.flag = TRUE: terminal-node
coefficients must decay and predictions must stay near the observed
response range. Fails against current source; the next commit fixes it.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: Apply the fix

**Files:**
- Modify: `boostmtree/R/boostmtree_fit_tree.R:741` and `:757`

**Interfaces:**
- Consumes: the failing test from Task 2.
- Produces: a converging in-sample boosting path under `cv.flag = TRUE`. No
  function signature changes; no exported behaviour renamed.

- [ ] **Step 1: Read the current control flow**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src/boostmtree
sed -n '714p;740,742p;755,758p' R/boostmtree_fit_tree.R
```

Expected exactly:

```
    if (model.info$cv.flag) {
      }
    } else {
      if (model.info$family == "nominal") {
        })
      }
    }
    for (q in seq_len(model.info$n.q)) {
```

If these lines differ, stop and re-derive the line numbers before editing.

- [ ] **Step 2: Turn the `else` into a plain close, and drop its trailing brace**

Line 741 becomes `    }` and line 757 is deleted:

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src/boostmtree
python3 - <<'PY'
p = "R/boostmtree_fit_tree.R"
lines = open(p).read().split("\n")

# 0-indexed: line 741 -> 740, line 757 -> 756
assert lines[740].strip() == "} else {", repr(lines[740])
assert lines[756].strip() == "}", repr(lines[756])

lines[740] = "    }"
del lines[756]

open(p, "w").write("\n".join(lines))
print("applied: else removed at 741, brace deleted at 757")
PY
```

- [ ] **Step 3: Verify the resulting structure and brace balance**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src/boostmtree
sed -n '712,716p;738,744p;753,758p' R/boostmtree_fit_tree.R
Rscript -e 'invisible(parse("R/boostmtree_fit_tree.R")); cat("parses cleanly\n")'
```

Expected: the CV block now closes with a bare `}`, the in-sample block
follows unconditionally, and `parses cleanly` is printed.

- [ ] **Step 4: Run the test and verify it PASSES**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R CMD INSTALL --no-docs --no-byte-compile -l "$BMLIB" boostmtree
R_LIBS_USER="$BMLIB" Rscript -e 'testthat::test_local("boostmtree", reporter="summary")'
```

Expected: PASS. `decay.ratio` ~0.204 and the predicted range ~`[-20.6, 4.5]`,
comfortably inside the allowed band.

- [ ] **Step 5: Commit**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
git add boostmtree/R/boostmtree_fit_tree.R
git commit -m "fix: refresh in-sample mu on every boosting iteration

The in-sample l.pred/mu update was the else branch of if (cv.flag), so
fitting with cv.flag = TRUE froze mu at its initialization. lambda and
rho/phi estimation both read that frozen mu, so the ridge penalty
collapsed and terminal-node coefficients never decayed, leaving the
stored gamma divergent and every non-CV prediction path unusable.

The CV state is additional bookkeeping; the in-sample linear predictor
must advance every iteration. Make the refresh unconditional.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Prove no collateral change to the `cv.flag = FALSE` path

**Files:**
- Create: `tools/baseline/compare.R`
- Create: `tools/baseline/patched-cvflag-false.rds` (generated)

**Interfaces:**
- Consumes: `tools/baseline/fingerprint.R` and
  `tools/baseline/baseline-cvflag-false.rds` from Task 1.
- Produces: `tools/baseline/compare.R`, runnable as
  `Rscript tools/baseline/compare.R <baseline.rds> <patched.rds>`; exits
  non-zero if any component differs.

- [ ] **Step 1: Write the comparison script**

Create `tools/baseline/compare.R`:

```r
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
```

- [ ] **Step 2: Generate the patched fingerprint**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R_LIBS_USER="$BMLIB" Rscript tools/baseline/fingerprint.R tools/baseline/patched-cvflag-false.rds
```

Expected: `fingerprint written to tools/baseline/patched-cvflag-false.rds`

- [ ] **Step 3: Compare and verify identical**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R_LIBS_USER="$BMLIB" Rscript tools/baseline/compare.R \
  tools/baseline/baseline-cvflag-false.rds \
  tools/baseline/patched-cvflag-false.rds
echo "exit=$?"
```

Expected: `cv.flag = FALSE path bit-identical across 5 components: mu, lambda, gamma, phi, rho` and `exit=0`.

If this reports a MISMATCH, the fix changed more than intended. Stop and
re-inspect the Task 3 diff with `git show HEAD -- boostmtree/R/boostmtree_fit_tree.R`.

- [ ] **Step 4: Commit**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
git add tools/baseline
git commit -m "test: verify cv.flag = FALSE path is unchanged by the fix

Bit-identical across mu, lambda, gamma, phi and rho.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Version bump, NEWS, and fork provenance

**Files:**
- Modify: `boostmtree/DESCRIPTION` (the `Version:`, `Date:` and `Description:` fields)
- Modify: `boostmtree/NEWS.md` (prepend a section)

**Note on line numbers:** Task 2 inserted two lines (`Suggests:` and
`Config/testthat/edition:`) after `Imports:`, so `Description:` is no longer
line 12. Locate these fields by name, not by line number.

**Interfaces:**
- Consumes: the fix from Task 3.
- Produces: `packageVersion("boostmtree") == "2.0.1"`.

- [ ] **Step 1: Bump version and date**

In `boostmtree/DESCRIPTION`, set the `Version:` field to `2.0.1` and the
`Date:` field to `2026-09-02`. These are lines 2 and 3 and were not moved by
Task 2, but match on the field name rather than the line number.

- [ ] **Step 2: Note the fork in the Description prose**

Append this sentence to the end of the `Description:` field, inside the same
field, after the existing final `<DOI:10.1007/s10994-016-5597-1>. `. Find the
field by matching `^Description:`, not by line number:

```
This is a Cleveland Clinic Foundation patched build correcting terminal-node coefficient divergence when fitting with cross-validation enabled.
```

- [ ] **Step 3: Verify DESCRIPTION still parses**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
Rscript -e 'd <- read.dcf("boostmtree/DESCRIPTION"); cat("Version:", d[,"Version"], "\n"); cat("fields:", ncol(d), "\n"); stopifnot(package_version(d[,"Version"]) == "2.0.1")'
```

Expected: `Version: 2.0.1` and no error.

- [ ] **Step 4: Prepend the NEWS entry**

At the very top of `boostmtree/NEWS.md`, above the existing `# boostmtree 2.0.0` line:

```markdown
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
  `predict(...)$muhat`, and `vimp()` on the non-CV path. Cross-validated
  output (`mu`, `err.rate`, `m.opt`, CV RMSE) was not affected.
* Added a regression test asserting the in-sample boost converges under
  `cv.flag = TRUE`.

```

- [ ] **Step 5: Commit**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
git add boostmtree/DESCRIPTION boostmtree/NEWS.md
git commit -m "chore: bump to 2.0.1 and document the cv.flag fix

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: `R CMD check` against the upstream baseline

**Files:**
- Create: `tools/baseline/patched-check.log` (generated)

**Interfaces:**
- Consumes: `tools/baseline/upstream-check.log` from Task 1.
- Produces: confirmation of no new WARNINGs or NOTEs.

- [ ] **Step 1: Run the check**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
R_LIBS_USER="$BMLIB" R CMD check --no-manual --no-build-vignettes boostmtree \
  > tools/baseline/patched-check.log 2>&1
tail -5 tools/baseline/patched-check.log
```

- [ ] **Step 2: Diff the problem lines against the baseline**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
diff <(grep -E "WARNING|NOTE|ERROR" tools/baseline/upstream-check.log | sort) \
     <(grep -E "WARNING|NOTE|ERROR" tools/baseline/patched-check.log | sort)
echo "exit=$?"
```

Expected: `exit=0`, or only additions you can explain and accept (a NOTE
about the newly added `tests/` directory is acceptable; a new WARNING is not).
If a new WARNING appears, fix it before proceeding.

- [ ] **Step 3: Confirm the test suite ran inside the check**

```bash
grep -A3 "checking tests" tools/baseline/patched-check.log
```

Expected: `OK`. If tests were skipped, testthat is not wired up — revisit Task 2.

- [ ] **Step 4: Clean up and commit the log**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
rm -rf boostmtree.Rcheck
git add tools/baseline/patched-check.log
git commit -m "test: record post-fix R CMD check log

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: Merge, tag, and publish install instructions

**Files:**
- Modify: none (git operations only)

**Interfaces:**
- Consumes: all prior tasks.
- Produces: tag `v2.0.1-ccf` on `ehrlinger/boostmtree_src`, installable by a third party.

- [ ] **Step 1: Push the branch and open the PR**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
git push -u origin fix/cv-flag-insample-mu
gh pr create --repo ehrlinger/boostmtree_src --base master \
  --title "fix: refresh in-sample mu on every boosting iteration (cv.flag = TRUE)" \
  --body "$(cat <<'EOF'
## Problem

Fitting with `cv.flag = TRUE` left the stored terminal-node coefficients
(`gamma`) divergent, so every output reading them returned values far
outside the observed response range: `predict()` with
`use.cv.flag = FALSE`, `partial.plot()`, `marginal.plot()`,
`predict(...)$muhat`, and `vimp()` on the non-CV path.

Measured on simulated data with response range `[-28.9, 6.4]`,
`predict()` returned `[-167, 82]` at `M = 200`. The divergence is linear
in `M`.

## Cause

`R/boostmtree_fit_tree.R:714` opened `if (model.info$cv.flag) {` and
`:741` made the in-sample branch its `else`, so the in-sample
`l.pred`/`mu` refresh ran only when `cv.flag = FALSE`. With CV enabled
`mu` stayed frozen at its initialization. `boostmtree.estimate.lambda()`
reads that frozen `mu` (under the default `control$cv.lambda = FALSE`), so
`lambda` collapsed by ~3 orders of magnitude, and `max|gamma|` never
decayed.

It went unnoticed because `finalize.fit` overwrites `fit$mu` with the
cross-validated mu, and the CV path is unaffected — CV error curves
differ by <0.1% between the broken and fixed builds.

## Fix

Make the in-sample refresh unconditional. The CV state is additional
bookkeeping; the in-sample linear predictor must advance every iteration.

## Verification

- `cv.flag = FALSE` results bit-identical to upstream across `mu`,
  `lambda`, `gamma`, `phi`, `rho`.
- New regression test asserts `max|gamma|` decays and predictions stay
  near the observed range under `cv.flag = TRUE`.
- `R CMD check` shows no new WARNINGs or NOTEs.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 2: Merge after review**

Merge the PR in the GitHub UI, or:

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
gh pr merge --repo ehrlinger/boostmtree_src --merge --delete-branch
git checkout master && git pull
```

- [ ] **Step 3: Tag the release**

```bash
cd /Users/ehrlinj/Documents/GitHub/boostmtree_src
git tag -a v2.0.1-ccf -m "boostmtree 2.0.1 (CCF patched build)

Fixes terminal-node coefficient divergence when fitting with cv.flag = TRUE."
git push origin v2.0.1-ccf
git rev-parse v2.0.1-ccf
```

Record the SHA — it is what goes in the manuscript methods statement.

- [ ] **Step 4: Verify a clean third-party install**

```bash
export CLEANLIB=/private/tmp/claude-504/-Users-ehrlinj-Documents-GitHub-boostmtree/5afbe653-7bf3-476e-b318-a4c4d253359b/scratchpad/cleanlib
mkdir -p "$CLEANLIB"
R_LIBS_USER="$CLEANLIB" Rscript -e '
remotes::install_github("ehrlinger/boostmtree_src", subdir = "boostmtree",
                        ref = "v2.0.1-ccf", lib = Sys.getenv("CLEANLIB"))
d <- packageDescription("boostmtree", lib.loc = Sys.getenv("CLEANLIB"))
cat("Version:", d$Version, "\n")
cat("RemoteUsername:", d$RemoteUsername, "\n")
cat("RemoteSha:", d$RemoteSha, "\n")'
```

Expected: `Version: 2.0.1`, `RemoteUsername: ehrlinger`, and a `RemoteSha`
matching the tag. This is the exact command the biostatistician runs.

- [ ] **Step 5: Refit the analysis models**

The corrupt `gamma` cannot be salvaged from existing `.cv` fit objects —
every EOA and mean-gradient model must be refit against 2.0.1 before any
figure or estimate goes into the manuscript. Both outcomes are separate
univariate fits (`n.q = 1`), so only the continuous path is involved.

---

## Post-plan follow-ons (not part of this plan)

- Upstream issue against `kogalur/boostmtree` with the minimal reproducer.
- hvtiPlotR helper for subgroup ensemble curves built from
  `predict(fit)$muhat` by row-subsetting, avoiding the `partial.plot()`
  `subset` -> `use.cv.flag = FALSE` downgrade entirely.
- `test-cv-flag-gamma-decay.R`'s regression coverage is `family = "continuous"`
  only. The gamma decay-ratio assertion as written is not portable to the
  binary, nominal, or ordinal families: their in-sample `mu` legitimately
  drives toward 0/1, so the ratio can stay above 1 even when the fix is
  correct. A family-appropriate assertion (or family-specific skip) is
  needed before extending this test to those families.
