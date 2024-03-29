% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rqf_funs.R
\name{rqfr}
\alias{rqfr}
\alias{rqfmr}
\alias{rqfp}
\title{Monte Carlo sampling of ratio/product of quadratic forms}
\usage{
rqfr(nit = 1000L, A, B, p = 1, q = p, mu, Sigma, use_cpp = TRUE)

rqfmr(nit = 1000L, A, B, D, p = 1, q = p/2, r = q, mu, Sigma, use_cpp = TRUE)

rqfp(nit = 1000L, A, B, D, p = 1, q = 1, r = 1, mu, Sigma, use_cpp = TRUE)
}
\arguments{
\item{nit}{Number of iteration or sample size.  Should be an integer-alike of
length 1.}

\item{A, B, D}{Argument matrices (see \dQuote{Details}).  Assumed to be square matrices of
the same order.  When missing, set to the identity matrix.  At least
one of these must be specified.}

\item{p, q, r}{Exponents for A, B, D, respectively (see \dQuote{Details}).  Assumed to be
numeric of length 1 each.  See \dQuote{Details} for default values.}

\item{mu}{Mean vector \eqn{\bm{\mu}}{\mu} for \eqn{\mathbf{x}}{x}.  Default
zero vector.}

\item{Sigma}{Covariance matrix \eqn{\mathbf{\Sigma}}{\Sigma} for
\eqn{\mathbf{x}}{x}.  Default identity matrix.  \code{mu} and
\code{Sigma} are assumed to be of the same order as the argument matrices.}

\item{use_cpp}{Logical to specify whether an \proglang{C++} version is called or
not.  \code{TRUE} by default.}
}
\value{
Numeric vector of length \code{nit}.
}
\description{
\code{rqfr()}, \code{rqfmr()}, and \code{rqfp()} calculate a random sample of
a simple ratio, multiple ratio (of special form), and product, respectively,
of quadratic forms in normal variables of specified mean and covariance
(standard multivariate normal by default).  These functions are primarily for
empirical verification of the analytic results provided in this package.
}
\details{
These functions generate a random sample of
\eqn{ \frac{(\mathbf{x^\mathit{T} A x})^p}{(\mathbf{x^\mathit{T} B x})^q}
     }{(x^T A x)^p / (x^T B x)^q}
(\code{rqfr()}),
\eqn{ \frac{(\mathbf{x^\mathit{T} A x})^p}
           {(\mathbf{x^\mathit{T} B x})^q (\mathbf{x^\mathit{T} Dx})^r}
     }{(x^T A x)^p / ( (x^T B x)^q (x^T D x)^r )}
(\code{rqfmr()}), and
\eqn{ (\mathbf{x^\mathit{T} A x})^p (\mathbf{x^\mathit{T} B x})^q
      (\mathbf{x^\mathit{T} D x})^r }{(x^T A x)^p (x^T B x)^q (x^T D x)^r}
(\code{rqfp()}), where
\eqn{\mathbf{x} \sim N_n(\bm{\mu}, \mathbf{\Sigma})
     }{x ~ N_n(\mu, \Sigma)}.  (Internally, \code{rqfr()} and \code{rqfmr()}
just call \code{rqfp()} with negative exponents.)

When only one of \code{p} and \code{q} is provided in \code{rqfr()},
the other (missing) one is set to the same value.

In \code{rqfmr()}, \code{q} and \code{r} are set to \code{p/2}
when both missing, and set to the same value when only one is missing.  When
\code{p} is missing, this is set to be \code{q + r}.  If unsure,
specify all these explicitly.

In \code{rqfp()}, \code{p}, \code{q} and \code{r} are \code{1} by default,
provided that the corresponding argument matrices are given.  If both
an argument matrix and its exponent (e.g., \code{D} and \code{r})
are missing, the exponent is set to \code{0} so that the factor be unity.
}
\examples{
p <- 4
A <- diag(1:p)
B <- diag(p:1)
D <- diag(sqrt(1:p))

## By default B = I, p = q = 1;
## i.e., (x^T A x) / (x^T x), x ~ N(0, I)
rqfr(5, A)

## (x^T A x) / ((x^T B x)(x^T D x))^(1/2), x ~ N(0, I)
rqfmr(5, A, B, D, 1, 1/2, 1/2)

## (x^T A x), x ~ N(0, I)
rqfp(5, A)

## (x^T A x) (x^T B x), x ~ N(0, I)
rqfp(5, A, B)

## (x^T A x) (x^T B x) (x^T D x), x ~ N(0, I)
rqfp(5, A, B, D)

## Example with non-standard normal
mu <- 1:p / p
Sigma <- matrix(0.5, p, p)
diag(Sigma) <- 1
rqfr(5, A, mu = 1:p / p, Sigma = Sigma)

## Compare Monte Carlo sample and analytic expression
set.seed(3)
mcres <- rqfr(1000, A, p = 2)
mean(mcres)
(anres <- qfrm(A, p = 2))
stats::t.test(mcres, mu = anres$statistic)

}
\seealso{
\code{\link{qfrm}} and \code{\link{qfpm}} for analytic moments

\code{\link{dqfr}} for analytic distribution-related functions for
simple ratios
}
