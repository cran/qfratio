% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dk_funs.R
\name{d1_i}
\alias{d1_i}
\alias{dtil1_i}
\alias{dtil1_i_v}
\alias{dtil1_i_m}
\title{Coefficients in polynomial expansion of generating function---single matrix}
\usage{
d1_i(L, m = 100L, thr_margin = 100)

dtil1_i_v(L, mu = rep.int(0, n), m = 100L, thr_margin = 100)

dtil1_i_m(A, mu = rep.int(0, n), m = 100L, thr_margin = 100)
}
\arguments{
\item{L}{Vector of eigenvalues of the argument matrix}

\item{m}{Integer-alike to specify the order of polynomials}

\item{thr_margin}{Optional argument to adjust the threshold for scaling
(see \dQuote{Details})}

\item{mu}{Mean vector \eqn{\bm{\mu}}{\mu} for \eqn{\mathbf{x}}{x}}

\item{A}{Argument matrix.  Assumed to be symmetric in these functions.}
}
\value{
Vector of length \code{m + 1}, corresponding to
the 0th, 1st, ..., and mth order terms.  Hence, the \code{[m + 1]}-th
element should be extracted when the coefficient for the \eqn{m}th order
term is required.

Has the attribute \code{"logscale"} as described in \dQuote{Scaling} above.
}
\description{
These are internal functions to calculate the coefficients
in polynomial expansion of generating functions for quadratic forms
in multivariate normal variables.

\code{d1_i()} is for standard multivariate normal variables,
\eqn{\mathbf{x} \sim N_n(\mathbf{0}_n, \mathbf{I}_n)}{x ~ N_n(0_n, I_n)}.

\code{dtil1_i_v()} is for noncentral multivariate normal variables,
\eqn{\mathbf{x} \sim N_n(\bm{\mu}, \mathbf{I}_n)}{x ~ N_n(\mu, I_n)}.

\code{dtil1_i_m()} is a wrapper for \code{dtil1_i_v()}
and takes the argument matrix rather than its eigenvalues.
}
\details{
\code{d1_i()} calculates \eqn{d_k(\mathbf{A})}{d_k(A)}, and
\code{dtil1_i_v()} and \code{dtil1_i_m()} calculate
\eqn{\tilde{d}_k(\mathbf{A})}{\tilde{d}_k(A)} in
Hillier et al. (2009, 2014) and Bao and Kan (2013).  The former is
related to the top-order zonal polynomial
\eqn{C_{[k]}(\mathbf{A})}{C_[k](A)} in the following way:
\eqn{ d_k(\mathbf{A}) = \frac{1}{k!} \left( \frac{1}{2} \right)_k
     C_{[k]}(\mathbf{A}) }{d_k(A) = (1/k!) (1/2)_k C_[k](A)},
where \eqn{(x)_k = x (x + 1) \dots (x + k - 1)}.

These functions calculate the coefficients based on the super-short
recursion algorithm described in Hillier et al. (2014: 3.2, eqs. 28--30).
\subsection{Scaling}{

The coefficients described herein (and in \code{\link{d2_ij}} and
\code{\link{d3_ijk}}) can become very large for higher-order terms,
so there is a practical risk of numerical overflow when applied to
large matrices or matrices with many large eigenvalues
(note that the latter typically arises from those with many small
eigenvalues for the front-end \code{qfrm()} functions).  To avoid
numerical overflow, these functions automatically scale
coefficients (and temporary objects used to calculate them) by a large number
(\code{1e10} at present) when any value in the temporary objects exceeds
a threshold, \code{.Machine$double.xmax / thr_margin / n}, where \code{n}
is the number of variables.  This default value empirically seems to work well
in most conditions, but use a large \code{thr_margin} (e.g., \code{1e5})
if you encounter numerical overflow.  (The \proglang{C++} functions use
an equivalent expression,
\code{std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n)}, with
\code{Scalar} being \code{double} or \code{long double}.)

In these \R functions, the scaling happens order-wise;
i.e., it influences all the coefficients of the same order in
multidimensional coefficients (in \code{\link{d2_ij}} and
\code{\link{d3_ijk}}) and the coefficients of the subsequent orders.

These scaling factors are recorded in the attribute \code{"logscale"} of the
return value, which is a vector/matrix/array whose size is identical to the
return value, so that \code{value / exp(attr(value, "logscale"))} equals
the original quantities to be obtained (if there were no overflow).

The \code{qfrm} and \code{qfmrm} functions handle return values of these
functions by first multiplying them with hypergeometric coefficients
(which are typically \eqn{\ll 1}{<< 1}) and then scaling the products back
to the original scale using the recorded scaling factors.  (To be precise,
this typically happens within \code{\link{hgs}} functions.)  The
\proglang{C++} functions handle the problem similarly (but by using
separate objects rather than attributes).

However, this procedure does not always mitigate the problem in
multiple series; when there are very large and very small
coefficients in the same order, smaller ones can diminish/underflow to
the numerical \code{0} after repeated scaling.  (The \code{qfrm} and
\code{qfmrm} functions try to detect and warn against
this by examining whether any of the highest-order terms is \code{0}.)
The present version of this package has implemented two methods to mitigate
this problem, but only through \proglang{C++} functions.  One is to use the
\code{long double} variable type, and the other is to use coefficient-wise
scaling (see \code{\link{qfrm}} and \code{\link{qfmrm}}).
}
}
\references{
Bao, Y. and Kan, R. (2013) On the moments of ratios of quadratic forms in
normal random variables. \emph{Journal of Multivariate Analysis}, \strong{117},
229--245.
\doi{10.1016/j.jmva.2013.03.002}.

Hillier, G., Kan, R. and Wang, X. (2009) Computationally efficient recursions
for top-order invariant polynomials with applications.
\emph{Econometric Theory}, \strong{25}, 211--242.
\doi{10.1017/S0266466608090075}.

Hillier, G., Kan, R. and Wang, X. (2014) Generating functions and
short recursions, with applications to the moments of quadratic forms
in noncentral normal vectors. \emph{Econometric Theory}, \strong{30}, 436--473.
\doi{10.1017/S0266466613000364}.
}
\seealso{
\code{\link{qfpm}}, \code{\link{qfrm}}, and \code{\link{qfmrm}} are
major front-end functions that utilize these functions

\code{\link{dtil2_pq}} for \eqn{\tilde{d}}
used for moments of a product of quadratic forms

\code{\link{d2_ij}} and \code{\link{d3_ijk}} for \eqn{d}, \eqn{h},
\eqn{\tilde{h}}, and \eqn{\hat{h}} used for moments of ratios
of quadratic forms
}
