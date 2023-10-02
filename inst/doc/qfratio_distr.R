## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, include = FALSE---------------------------------------------------
library(qfratio)
require(stats)
set.seed(764561)

## ----definition_exported, eval = FALSE----------------------------------------
#  pqfr <- function(quantile, A, B, p = 1, mu = rep.int(0, n), Sigma = diag(n),
#                   lower.tail = TRUE, log.p = FALSE,
#                   method = c("imhof", "davies", "forchini", "butler"),
#                   trim_values = TRUE, return_abserr_attr = FALSE, m = 100L,
#                   tol_zero = .Machine$double.eps * 100,
#                   tol_sing = tol_zero, ...) { ... }
#  dqfr <- function(quantile, A, B, p = 1, mu = rep.int(0, n), Sigma = diag(n),
#                   log = FALSE, method = c("broda", "hillier", "butler"),
#                   trim_values = TRUE, normalize_spa = FALSE,
#                   return_abserr_attr = FALSE, m = 100L,
#                   tol_zero = .Machine$double.eps * 100,
#                   tol_sing = tol_zero, ...) { ... }

## ----example_methods, error = TRUE--------------------------------------------
## Choice from alternative methods
A <- diag(1:3)
pqfr(1.5, A, method = "imhof")    # default
pqfr(1.5, A, method = "davies")   # similar
pqfr(1.5, A, method = "forchini") # series
pqfr(1.5, A, method = "butler")   # spa

dqfr(1.5, A, method = "broda")    # default
dqfr(1.5, A, method = "hillier")  # series
dqfr(1.5, A, method = "butler")   # spa

## Not recommended; for diagnostic use only
qfratio:::pqfr_imhof(1.5, A)
qfratio:::pqfr_A1B1(1.5, A, m = 9, check_convergence = FALSE)


## This is okay
x <- c(1.5, 2.5, 3.5)
pqfr(x, A)

## This is not
qfratio:::pqfr_imhof(x, A)

## ----ks_syntax_correct--------------------------------------------------------
## Small Monte Carlo sample
A <- diag(1:3)
B <- diag(sqrt(1:3))
x <- rqfr(10, A, B)

## Calculate p-values
pseq <- pqfr(x, A, B, return_abserr_attr = TRUE)

## Maximum error when evaluated at x;
## looks small enough
max(attr(pseq, "abserr"))

## Correct syntax, expected outcome
## \(q) syntax could also be used in recent versions of R
ks.test(x, function(q) pqfr(q, A, B))

## ----ks_syntax_wrong----------------------------------------------------------
## Incorrect; no error/warning because
## B is passed to ks.test rather than to pqfr
ks.test(x, pqfr, A = A, B = B)

## ----definition_internal_series, eval = FALSE---------------------------------
#  ## Used in pqfr(..., method = "forchini")
#  pqfr_A1B1 <- function(quantile, A, B, m = 100L, mu = rep.int(0, n),
#                        check_convergence = c("relative", "strict_relative",
#                                              "absolute", "none"),
#                        stop_on_error = FALSE, use_cpp = TRUE,
#                        cpp_method = c("double", "long_double", "coef_wise"),
#                        nthreads = 1,
#                        tol_conv = .Machine$double.eps ^ (1/4),
#                        tol_zero = .Machine$double.eps * 100,
#                        tol_sing = tol_zero,
#                        thr_margin = 100) { ... }
#  ## Used in dqfr(..., method = "hillier")
#  dqfr_A1I1 <- function(quantile, LA, m = 100L,
#                        check_convergence = c("relative", "strict_relative",
#                                              "absolute", "none"),
#                        use_cpp = TRUE,
#                        tol_conv = .Machine$double.eps ^ (1/4),
#                        thr_margin = 100) { ... }

## ----example_series_1, error = TRUE-------------------------------------------
A <- diag(1:3)
pqfr(1.5, A, method = "forchini")
dqfr(1.5, A, method = "hillier")

B <- diag(sqrt(1:3))
pqfr(1.5, A, B, method = "forchini")
## dqfr method does not accommodate B, mu, or Sigma
dqfr(1.5, A, B, method = "hillier")

## ----example_series_2, error = TRUE-------------------------------------------
A <- diag(1:3)

## p-value just below 2, an eigenvalue of A
## Typically throws two warnings:
##   Maximum iteration in hypergeometric function
##   and non-convergence of series
pqfr(1.9999, A, method = "forchini")

## More realistic value; expected from symmetry
pqfr(1.9999, A, method = "imhof")

## ----definition_internal_inversion, eval = FALSE------------------------------
#  ## Used in pqfr(..., method = "imhof") (default)
#  pqfr_imhof <- function(quantile, A, B, mu = rep.int(0, n),
#                         autoscale_args = 1, stop_on_error = TRUE, use_cpp = TRUE,
#                         tol_zero = .Machine$double.eps * 100,
#                         epsabs = epsrel, epsrel = 1e-6, limit = 1e4) { ... }
#  ## Used in pqfr(..., method = "davies")
#  pqfr_davies <- function(quantile, A, B, mu = rep.int(0, n),
#                          autoscale_args = 1,
#                          tol_zero = .Machine$double.eps * 100, ...) { ... }
#  ## Used in dqfr(..., method = "broda") (default)
#  dqfr_broda <- function(quantile, A, B, mu = rep.int(0, n),
#                         autoscale_args = 1, stop_on_error = TRUE,
#                         use_cpp = TRUE, tol_zero = .Machine$double.eps * 100,
#                         epsabs = epsrel, epsrel = 1e-6, limit = 1e4) { ... }

## ----example_specify_error----------------------------------------------------
A <- diag(1:4)

## This error bound satisfies "abserr < value * epsrel"
pqfr(3.9, A, method = "imhof", return_abserr_attr = TRUE,
     epsabs = 0, epsrel = 1e-6)

## This one violates "abserr < value * epsrel",
## although abserr is a valid error bound
pqfr(1.2, A, method = "imhof", return_abserr_attr = TRUE,
     epsabs = 0, epsrel = 1e-6)


## ----example_inversion_scale, error = TRUE------------------------------------
A <- diag(1:3)
B <- diag(sqrt(1:3))

## Without autoscale_args
## We know these are equal
pqfr(1.5, A, B, autoscale_args = FALSE)
pqfr(1.5, A * 1e-10, B * 1e-10, autoscale_args = FALSE)
## The latter failed because of numerically small eigenvalues

## With autoscale_args = 1 (default)
pqfr(1.5, A * 1e-10, B * 1e-10)

## ----example_inversion_trim, error = TRUE-------------------------------------
## Result without trimming;
## (typically) negative density, which is absurd
## In this case, error interval typically spans across 0
dqfr(1.2, diag(1:30), return_abserr_attr = TRUE,
     trim_values = FALSE)

## Result with trimming (default)
dqfr(1.2, diag(1:30), return_abserr_attr = TRUE)
## Note that the actual value is only bounded by
## 0 and abserr

## ----definition_internal_spa, eval = FALSE------------------------------------
#  ## Used in pqfr(..., method = "butler")
#  pqfr_butler <- function(quantile, A, B, mu = rep.int(0, n),
#                          order_spa = 2, stop_on_error = FALSE, use_cpp = TRUE,
#                          tol_zero = .Machine$double.eps * 100,
#                          epsabs = .Machine$double.eps ^ (1/2), epsrel = 0,
#                          maxiter = 5000) { ... }
#  ## Used in dqfr(..., method = "butler")
#  dqfr_butler <- function(quantile, A, B, mu = rep.int(0, n),
#                          order_spa = 2, stop_on_error = FALSE, use_cpp = TRUE,
#                          tol_zero = .Machine$double.eps * 100,
#                          epsabs = .Machine$double.eps ^ (1/2), epsrel = 0,
#                          maxiter = 5000) { ... }

## ----example_spa_1------------------------------------------------------------
A <- diag(1:3)

## Default for spa distribution function
pqfr(1.2, A, method = "butler", order_spa = 2)

## First-order spa
pqfr(1.2, A, method = "butler", order_spa = 1)

## More accurate numerical inversion
pqfr(1.2, A)


## Default for density
dqfr(1.2, A, method = "butler",
     order_spa = 2, normalize_spa = FALSE)

## First-order
dqfr(1.2, A, method = "butler",
     order_spa = 1, normalize_spa = FALSE)

## Normalized density, second-order
dqfr(1.2, A, method = "butler",
     order_spa = 2, normalize_spa = TRUE)

## Normalized density, first-order
dqfr(1.2, A, method = "butler",
     order_spa = 1, normalize_spa = TRUE)

## More accurate numerical inversion
dqfr(1.2, A)


## ----errorbound_1-------------------------------------------------------------
A <- diag(1:4)

pqfr(1.5, A, return_abserr_attr = TRUE)

dqfr(1.5, A, return_abserr_attr = TRUE)

## ----errorbound_trim, error = TRUE--------------------------------------------
## Without trimming, result is (typically) negative
## But note that value + abserr is positive
dqfr(1.2, diag(1:35), return_abserr_attr = TRUE,
     epsabs = 1e-10, trim_values = FALSE)

## With trimming, value is replaced by tol_zero
## Note slightly shortened abserr
dqfr(1.2, diag(1:35), return_abserr_attr = TRUE,
     epsabs = 1e-10)


## When untrimmed value + abserr < tol_zero
dqfr(1.1, diag(1:35), return_abserr_attr = TRUE,
     epsabs = 1e-15, trim_values = FALSE)
## True value is somewhere between 0 and value + abserr
## (assuming these are reliable)

## When trimmed, abserr reflects tol_zero
## because the true value is between 0 and tol_zero
dqfr(1.1, diag(1:35), return_abserr_attr = TRUE,
     epsabs = 1e-15)

## ----example_profile_distr, fig.width = 4, figh.height = 4, error = TRUE------
A <- diag(1:4)
qseq <- seq.int(0.8, 4.2, length.out = 100)

## Generate p-value sequences
## Warning is expected
pseq_inv <- pqfr(qseq, A, method = "imhof",
                 return_abserr_attr = TRUE)
pseq_ser <- pqfr(qseq, A, method = "forchini",
                 check_convergence = FALSE)
pseq_spa <- pqfr(qseq, A, method = "butler")

## Maximum error in numerical inversion;
## looks small enough
max(attr(pseq_inv, "abserr"))

## Graphical comparison
par(mar = c(4, 4, 0.1, 0.1))
plot(qseq, type = "n", xlim = c(1, 4), ylim = c(0, 1),
     xlab = "q", ylab = "F(q)")
lines(qseq, pseq_inv, col = "gray", lty = 1)
lines(qseq, pseq_ser, col = "tomato", lty = 2)
lines(qseq, pseq_spa, col = "slateblue", lty = 3)
legend("topleft", legend = c("inversion", "series", "saddlepoint"),
       col = c("gray", "tomato", "slateblue"), lty = 1:3, cex = 0.8)

## Logical vector to exclude q around eigenvalues of A
avoid_evals <- ((qseq %% 1) > 0.05) & ((qseq %% 1) < 0.95)

## Numerical comparison
all.equal(pseq_inv[avoid_evals], pseq_ser[avoid_evals],
          check.attributes = FALSE)
all.equal(pseq_inv[avoid_evals], pseq_spa[avoid_evals],
          check.attributes = FALSE)

## ----example_profile_density, fig.width = 4, figh.height = 4, error = TRUE----
## Generate p-value sequences
dseq_inv <- dqfr(qseq, A, method = "broda",
                 return_abserr_attr = TRUE)
dseq_ser <- dqfr(qseq, A, method = "hillier",
                 check_convergence = FALSE)
dseq_spa <- dqfr(qseq, A, method = "butler")

## Maximum error in numerical inversion;
## looks small enough
max(attr(dseq_inv, "abserr"))

## Graphical comparison
par(mar = c(4, 4, 0.1, 0.1))
plot(qseq, type = "n", xlim = c(1, 4), ylim = c(0, 0.8),
     xlab = "q", ylab = "f(q)")
lines(qseq, dseq_inv, col = "gray", lty = 1)
lines(qseq, dseq_ser, col = "tomato", lty = 2)
lines(qseq, dseq_spa, col = "slateblue", lty = 3)
legend("topleft", legend = c("inversion", "series", "saddlepoint"),
       col = c("gray", "tomato", "slateblue"), lty = 1:3, cex = 0.8)

## Numerical comparison
all.equal(dseq_inv, dseq_ser, check.attributes = FALSE)
all.equal(dseq_inv, dseq_spa, check.attributes = FALSE)

## Do densities sum up to 1?
sum(dseq_inv * diff(qseq)[1])
sum(dseq_ser * diff(qseq)[1])
sum(dseq_spa * diff(qseq)[1])

## ----example_density_normalize, fig.width = 4, figh.height = 4, error = TRUE----
## Normalized saddlepoint approximation density
dseq_spa_normalized <- dqfr(qseq, A, method = "butler",
                            normalize_spa = TRUE)
all.equal(dseq_inv, dseq_spa_normalized,
          check.attributes = FALSE)
sum(dseq_spa_normalized * diff(qseq)[1])

