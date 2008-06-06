\name{gam.method}
\alias{gam.method}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Setting GAM fitting method}
\description{ This is a function of package \code{mgcv} which allows 
selection of the numerical method used to optimize the  smoothing parameter 
estimation criterion for a gam. It is used to set argument \code{method} of \code{\link{gam}}.
}
\usage{
gam.method(gam="outer",outer="newton",gcv="deviance",
           family=NULL)
}
\arguments{ 

\item{gam}{Which method to use in the generalized case (i.e. all case other
than gaussian with identity link). \code{"perf"} for the performance
iteration (see details) with \code{\link{magic}} as the basic estimation
engine. \code{"perf.outer"} for \code{\link{magic}} based performance
iteration followed by outer iteration (see details). \code{"outer"} for pure
outer iteration.}

\item{outer}{ The optimization approach to use to optimize log smoothing
parameters by outer
iteration.\code{"newton"} (default) for modified Newton method backed up by steepest
descent, based on exact first and second derivatives. \code{"bfgs"} for a hybrid Newton-
Quasi-Newton approach which can be faster for models with many smoothing parameters.
 \code{"nlm"} to use 
\code{\link{nlm}} with exact first derivatives
to optimize the smoothness selection criterion. \code{"nlm.fd"} to use
\code{nlm} with finite differenced first derivatives (slower and less
reliable). \code{"optim"} to use the \code{"L-BFGS-B"} quasi-Newton method
option of routine \code{optim}, with exact first derivatives.}

\item{gcv}{One of \code{"deviance"} or \code{"GACV"},
specifying the flavour of GCV to use with outer iteration. \code{"deviance"}
simply replaces the residual sum of squares term in a GCV score with the
deviance, following Hastie and Tibshirani (1990, section 6.9). \code{"GACV"}
(only available with outer method \code{"newton"}) 
uses a varient of Gu and Xiang's (2001) generalized approximate cross
validation, modified to deal with arbitrary link-error
combinations (see Wood, 2008). 
}

\item{family}{The routine is called by \code{\link{gam}} to check the supplied
method argument. In this circumstance the family argument is passed, to check
that it works with the specified method. }
}

\details{ The default methods used by \code{gam} are based on Newton type optimization of
GCV/UBRE/AIC scores with respect to smoothing parameters, as described in Wood (2004) 
and Wood (2008) for the additive and generalized additive model cases respectively. 
The smoothing criteria (GCV etc) are evaluated for the fitted model itself (rather than some working 
approximate model). Since this involves optimizing the criteria `outside' the PIRLS or LS method used 
for fitting, it is referred to as `outer' iteration. 

In the generalized case several alternative optimisation methods can be used for outer
optimization. Usually the fastest and most
reliable approach is to use a modified Newton optimizer with exact first and
second derivatives, and this is the default. However if there are large numbers of smoothing parameters
then it can be faster to use a hybrid Newton-BFGS quasi-Newton approach, which avoids the expense of 
frequent full second derivative evaluation.

 \code{nlm} can be used with finite differenced first derivatives. This is not ideal theoretically, since 
it is possible for the finite difference estimates of derivatives to be very
badly in error on rare
occasions when the P-IRLS convergence tolerance is close to being matched
exactly, so that two components of a finite differenced derivative require
different numbers of iterations of P-IRLS in their evaluation. An alternative
is provided in which \code{nlm} uses numerically exact first derivatives, this
is faster and less problematic than the other scheme. A further alternative is to use a quasi-Newton
scheme with exact derivtives, based on \code{optim}. In practice this usually
seems to be slower than the \code{nlm} method. 


The alternative approach of `performance oriented iteration' was suggested by Gu 
(and is rather similar to the PQL method in generalized linear mixed modelling). At each step of the
P-IRLS (penalized iteratively reweighted least squares) iteration, 
by which a gam is fitted, the smoothing parameters are
estimated by GCV or UBRE applied to the working penalized linear modelling
problem. In most cases, this process converges and gives smoothness estimates
that perform well (see e.g. Wood 2004). 

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

In summary: performance iteration can fail to converge. It may occasionally be 
faster than outer iteration, but since the Wood (2008) method `outer' iteration 
is faster in most cases.
}


\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass.

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for generalized
additive models. J.R.Statist.Soc.B 70(3):495-518

\url{http://www.maths.bath.ac.uk/~sw283/}


}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{\code{\link{gam.control}}   \code{\link{gam}}, \code{\link{gam.fit}}, \code{\link{glm.control}} }

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..




