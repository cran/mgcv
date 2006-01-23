\name{gam.method}
\alias{gam.method}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Setting GAM fitting method}
\description{ This is a function of package \code{mgcv} which allows 
selection of the numerical method used to optimize the  smoothing parameter 
estimation criterion for a gam.

It is used to set argument \code{method} of \code{\link{gam}}.
}
\usage{
gam.method(am="magic",gam="outer",outer="nlm",pearson=FALSE,family=NULL)
}
\arguments{ 
\item{am}{ Which method to use for a pure additive model (i.e. identity link,
gaussian errors). Either \code{"magic"} if the Wood (2004) method
(\code{\link{magic}}) is to be used, or \code{"mgcv"} if the faster, but less
stable and general Wood (2000) method (\code{\link{mgcv}}) should be used.}

\item{gam}{Which method to use in the generalized case (i.e. all case other
than gaussian with identity link). \code{"perf.magic"} for the performance
iteration (see details) with \code{\link{magic}} as the basic estimation
engine. \code{"perf.mgcv"} for performance iteration with
\code{\link{mgcv}} as the underlying estimation
engine. \code{"perf.outer"} for \code{\link{magic}} based performance
iteration followed by outer iteration (see details). \code{"outer"} for pure
outer iteration.}

\item{outer}{ The optimization approach to use for outer
iteration. \code{"nlm"} to use \code{\link{nlm}} with exact first derivatives
to optimize the smoothness selection criterion. \code{"nlm.fd"} to use
\code{nlm} with finite differenced first dericatives (slower and less
reliable). \code{"optim"} to use the \code{"L-BFGS-B"} quasi-Newton method
option of routine \code{optim}, with exact first derivatives.}

\item{pearson}{When using any sort of outer iteration, there is a choise
between using GCV/UBRE scores based on the Pearson statistic, or based on the
deviance. The Pearson versions result in smoother models, but this is often
oversmoothing, particularly for `low count' data and binary regression.}

\item{family}{The routine is called by \code{\link{gam}} to check the supplied
method argument. In this circumstance the family argument is passed, to check
that it works with the specified method. Negative binomial families only work 
with performance iteration, and the method is reset to this if necessary.}
}

\details{ The performance iteration was suggested by Gu (and is rather similar
to the PQL method in generalized linear mixed modelling). At each step of the
P-IRLS (penalized iteratively reweighted least squares) iteration, 
by which a gam is fitted, the smoothing parameters are
estimated by GCV or UBRE applied to the working penalized linear modelling
problem. In most cases, this process converges and gives smoothness estimates
that perform well. It is usually very fast, since smoothing parameters are
estimated alongside other model coefficients in a single P-IRLS process.

The performance iteration has two disadvantages. (i) in the presence of
co-linearity or concurvity (a frequent problem when spatial smoothers are
included in a model with other covariates) then the process can fail to
converge. Suppose we start with some coefficient and smoothing parameter
estimates, implying a working penalized linear model: the optimal smoothing
parameters and coefficients for this working model may in turn imply a working 
model for which the original estimates are better than the most recent
estimates. This sort of effect can prevent convergence.

Secondly it is often possible to find a set of smoothing parameters that
result in a lower GCV or UBRE score, for the final working model, than the
final score that results from the performance iterations. This is because the
performance iteration is only approximately optimizing this score (since
optimization is only performed on the working model). The disadvantage here is
not that the model with lower score would perform better (it usually doesn't),
but rather that it makes model comparison on the basis of GCV/UBRE score
rather difficult. 

Both disadvantages of performance iteration are surmountable by using what is
basically O'Sullivan's (1986) suggestion. Here the P-IRLS scheme is iterated
to convergence for a fixed set of smoothing parameters, with an appropriate
GCV/UBRE score evaluated at convergence. This score at convergence is
optimized in some way. This is termed \code{"outer"} optimization, since the
optimization is outer to the P-IRLS loop. Outer iteration is slower than
performance iteration.

The `appropriate GCV/UBRE' score in the previous paragraph can be defined in
one of two ways either (i) the deviance, or (ii) the Pearson statistic can be
used in place of the residual sum of squares in the GCV/UBRE score. (ii) makes
the GCV/UBRE score correspond to the score for the working linear model at
convergence of the P-IRLS, but in practice tends to result in oversmoothing, 
particularly with low n binomial data, or low mean counts. Hence the default
is to use (i).

Several alternative optimisation methods can be used for outer
optimization. \code{nlm} can be used with finite differenced first
derivatives. This is not ideal theoretically, since it is possible for the
finite difference estimates of derivatives to be very badly in error on rare
occasions when the P-IRLS convergence tolerance is close to being matched
exactly, so that two components of a finite differenced derivative require
different numbers of iterations of P-IRLS in their evaluation. An alternative
is provided in which \code{nlm} uses numerically exact first derivatives, this
is faster and less problematic than the other scheme. Finally, a quasi-Newton
scheme with exact derivtives can be used instead, based on \code{optim}. In practice this usually
seems to be slower than the \code{nlm} method.  

It is possible to iterate the performance iteration to convergence and then
improve the smoothing parameter estimates using outer iteration: only a few steps are
usually required in the outer iteration in this case, so it may be quite
efficient, but it is not recommended if the performance iteration itself is
non-convergent. When using `pure' outer iteration, a single step of the
performance iteration is in fact taken first, to obtain estimates of the scale
of the GCV/UBRE objective: starting values for the smoothing parameters are
obtained using \code{\link{initial.sp}}.

In summary: performance iteration is fast, but can fail to converge. Outer
iteration is slower, but more reliable. At present only performance iteration
is available for negative binomial families.
}


\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 
65(1):95-114

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass.


\url{http://www.maths.bath.ac.uk/~sw283/}


}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{\code{\link{gam.control}}   \code{\link{gam}}, \code{\link{gam.fit}}, \code{\link{glm.control}} }

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..




