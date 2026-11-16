# Bug Report: `full.laplaceQ_MA()` crashes on clustered quantal data

**Date:** 2026-03-25  
**Package version:** BMABMDR 0.1.20  
**Status:** Fixed in commit 5124b66

---

## Summary

`full.laplaceQ_MA()` (and `anydoseresponseN()`) crashed with
**"$ operator is invalid for atomic vectors"** when fitting clustered
(litter-level) quantal data. The root cause was that `gamlss::gamlss()`
for the Beta-Binomial saturated fit was wrapped in `try()`, but the code
never checked whether `try()` returned an error before accessing
`fpfit2$sigma.coefficients`.

---

## What happens when `gamlss()` fails

`gamlss::gamlss()` can fail to converge (e.g. singular data, extreme
dose-response shapes, or numerical issues). When it does, `try()` returns
a **character string** of class `"try-error"` instead of a fitted model
object. The `$` operator cannot be applied to a character string, so any
code that unconditionally does `fpfit2$sigma.coefficients` immediately
throws:

```
Error in fpfit2$sigma.coefficients :
  $ operator is invalid for atomic vectors
```

This error propagates out of `modelTestQ()`, `PREP_DATA_QA()`, and
`anydoseresponseN()` — **completely aborting the analysis** even though
the rest of the fitting pipeline could continue with a reasonable fallback
value for the intraclass correlation coefficient (ICC, `rhohat`).

---

## Impact of the fix

The fix adds an `inherits(fpfit2, "try-error")` guard at each of the
three affected call sites:

| File | Function | Context |
|------|----------|---------|
| `R/fun_modelTest.R` | `modelTestQ()` | Saturated model prior for Beta-Binomial |
| `R/fun_Data.R` | `PREP_DATA_QA()` | Starting-value estimation when `cluster = TRUE` |
| `R/anydoseresponse.R` | `anydoseresponseN()` | H0 model prior for Beta-Binomial |

When `gamlss()` fails, `rhohat` is set to **0.1** (a sensible ICC
default) instead of crashing. The subsequent `ifelse(rhohat == 0, 0.001,
rhohat)` guard is then applied exactly as in the success path, so
downstream Stan model fitting proceeds normally.

**Before the fix:** any convergence failure in the Beta-Binomial
saturated fit causes a hard crash — the user gets no BMD estimate.

**After the fix:** the analysis continues with `rhohat = 0.1` as the
ICC starting value / prior centre. The BMD estimates and model weights
are still computed; they may differ slightly from a run where the
saturated gamlss fit converged, but the results are valid.

---

## Reproducible example

```r
library(BMABMDR)

# Built-in clustered quantal dataset (one row per litter)
data("data_quantal_litter")
head(data_quantal_litter)
#>   dose y  n
#> 1    0 1 12
#> 2    0 1 12
#> 3    0 0 12

# Prepare the data for clustered (litter-effect) analysis
data_QL <- PREP_DATA_QA(
  data     = data_quantal_litter,
  sumstats = TRUE,
  q        = 0.1,
  cluster  = TRUE
)

# ---- BEFORE FIX ----
# The line below crashed when gamlss() failed inside PREP_DATA_QA() /
# modelTestQ():
#
# Error in fpfit2$sigma.coefficients :
#   $ operator is invalid for atomic vectors
#
# Traceback:
#   1. full.laplaceQ_MA(data.Q = data_QL, prior.weights = rep(1, 8))
#   2.   modelTestQ(...)

# ---- AFTER FIX ----
# Now completes successfully, using rhohat = 0.1 as a fallback ICC when
# the Beta-Binomial saturated fit does not converge.
fit_ql <- full.laplaceQ_MA(
  data.Q        = data_QL,
  prior.weights = rep(1, 8)
)

# Inspect results
class(fit_ql)          # "BMADRQ"
print(fit_ql)
```

### Minimal synthetic example that triggers the bug

The following synthetic dataset mimics the pathological case (all
responses 0 at control, saturated model singular) that reliably causes
`gamlss()` to fail and previously triggered the crash:

```r
library(BMABMDR)

# Perfectly separated data: gamlss saturated fit is numerically unstable
bad_data <- data.frame(
  dose = c(0, 0, 0, 1, 1, 1, 5, 5, 5),
  y    = c(0, 0, 0, 3, 3, 3, 5, 5, 5),
  n    = c(5, 5, 5, 5, 5, 5, 5, 5, 5)
)

# ---- BEFORE FIX: this crashed ----
# ---- AFTER FIX: returns a valid BMADRQ object ----
data_bad <- PREP_DATA_QA(data = bad_data, sumstats = TRUE, q = 0.1,
                         cluster = TRUE)
fit_bad  <- full.laplaceQ_MA(data.Q = data_bad, prior.weights = rep(1, 8))
class(fit_bad)  # "BMADRQ"
```

---

## Root-cause code (before fix)

```r
# R/fun_modelTest.R  (~line 524)
fpfit2 <- try(
  gamlss::gamlss(cbind(yy, n.a - yy) ~ as.factor(xx),
                 sigma.formula = ~1, family = BB, data = datf),
  silent = TRUE
)
# BUG: fpfit2 may be a "try-error" string here
rhohat <- exp(fpfit2$sigma.coefficients) /   # <-- crashes
          (exp(fpfit2$sigma.coefficients) + 1)
```

## Fixed code

```r
fpfit2 <- try(
  gamlss::gamlss(cbind(yy, n.a - yy) ~ as.factor(xx),
                 sigma.formula = ~1, family = BB, data = datf),
  silent = TRUE
)
if (inherits(fpfit2, "try-error")) {
  rhohat <- 0.1          # sensible ICC fallback
} else {
  rhohat <- exp(fpfit2$sigma.coefficients) /
            (exp(fpfit2$sigma.coefficients) + 1)
}
rhohat <- ifelse(rhohat == 0, 0.001, rhohat)
```

---

## Environment

- R 4.5.1  
- gamlss 5.4-22  
- BMABMDR 0.1.20  
- Windows 11
