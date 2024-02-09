## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, include = FALSE---------------------------------------------------
library(qfratio)
set.seed(764561)

## ----first_example------------------------------------------------------------
## Simple matrices
nv <- 4
A <- diag(1:nv)
B <- diag(sqrt(nv:1))
D <- diag((nv:1) ^ 2)

## Expectation of (x^T A x)^2 / (x^T x)^1 where x ~ N(0, I)
qfrm(A, p = 2, q = 1) # By default, B = I

## Expectation of (x^T A x)^2 / (x^T B x)^2 where x ~ N(0, I)
qfrm(A, B, p = 2) # By default, q = p

## Expectation of (x^T A x)^1/2 / (x^T B x)^1 where x ~ N(0, I)
qfrm(A, B, p = 1/2, q = 1)

## Expectation of (x^T A x)^2 / ((x^T B x)^1 (x^T D x)^1) where x ~ N(0, I)
qfmrm(A, B, D, p = 2, q = 1, r = 1)

## Expectation of (x^T A x)^1 / (x^T B x)^1 where x ~ N(mu, Sigma)
## with arbitrary mu and Sigma
mu <- 1:nv / nv
Sigma <- diag(runif(nv) * 3)
qfrm(A, B, p = 1, q = 1, mu = mu, Sigma = Sigma, m = 300)

## ----example, fig.width = 4, fig.height = 4-----------------------------------
## Preparing simple objects
nv <- 10
A <- diag(1:nv)
B <- diag(sqrt(nv:1))

#####
## I: Simple case, (x^T A t)^2 / (x^T x)^2, where x ~ N(0, I)
##    Exact solution available
res <- qfrm(A, p = 2, q = 2)
res
## In this case, plot is moot because the result is exact

#####
## II: Simple case, (x^T A t)^2 / (x^T B x)^2, where x ~ N(0, I)
##     Error bound available
res <- qfrm(A, B, 2, 2)
res

## Displaying more digits
print(res, digits = 10)

## Direct access to the results:
names(res)
res$statistic       # Evaluation result (partial sum)
res$error_bound     # Error bound for $statistic
tail(res$terms)     # Tail of truncated series; sum($terms) == $statistic
tail(res$seq_error) # Tail of sequence of error bounds; last one is $error_bound

## Inspect plot
plot(res)
## Note that profile for error bound is shown as
## upper bound of the true value in the plot

#####
## III: Simple case, (x^T A t)^(1/2) / (x^T x)^(1/2), where x ~ N(0, I)
##      Error bound can be negative
res <- qfrm(A, p = 1/2, q = 1/2)
res # Note negative (though one-sided) error bound

## Inspect plot
plot(res)

#####
## IV: Non-central case, (x^T A t)^2 / (x^T B x)^2, where x ~ N(mu, I)
##     Two-sided error bound
mu <- 1:nv / sum(1:nv)
res <- qfrm(A, B, p = 2, q = 2, mu = mu)
res # Note two-sided error bound

## plot() automatically handles two-sided error bound
plot(res)


#####
## V: Difficult case with highly dispersed eigenvalues
##    Error bound available
A <- diag((1:nv) ^ 2)
B <- diag((nv:1) ^ 3)

res <- qfrm(A, B, 2, 2) # Note the warning
res

plot(res) # Non-convergence evident

## Try larger m
res <- qfrm(A, B, 2, 2, m = 2000)
res

plot(res) # Much better


#####
## VI: Difficult case with highly dispersed eigenvalues
##     Error bound unavailable

## It is recommended to use a tol_conv stricter than
## the default value, which is .Machine$double.eps ^ (1/4)
res <- qfrm(A, B, 1/2, 1/2, tol_conv = sqrt(.Machine$double.eps))
res

plot(res) # Close to convergence, but still growing

## Try larger m
res <- qfrm(A, B, 1/2, 1/2, m = 500, tol_conv = sqrt(.Machine$double.eps))
res

plot(res) # Better

## ----MC, fig.width = 4, fig.height = 4----------------------------------------
## A large problem
nv <- 200
large_A <- diag(c(1000, rep.int(1, nv - 1)))
large_B <- diag(c(rep.int(1, nv - 1), 1000))
large_D <- diag((nv:1) ^ 2)

res <- qfmrm(large_A, large_B, large_D, 1, 1/2, 1/2, m = 500)
res

plot(res) # Far to convergence

## This problem needs much larger m (>5000)
## taking computational time + memory
## (though still manageable with a regular machine in this particular case)

## Monte Carlo sample
MCres <- rqfmr(10000, large_A, large_B, large_D, 1, 1/2, 1/2)

## Monte Carlo 95% CI of the moment
mean(MCres) + sd(MCres) / sqrt(10000) * qt(c(0.025, 0.975), 10000 - 1)

