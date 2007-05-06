\name{gam.fit2}
\alias{gam.fit2}
\alias{gam.fit3}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{P-IRLS GAM estimation with GCV \& UBRE derivative calculation}
\description{Estimation of GAM smoothing parameters is most stable if
optimization of the UBRE or GCV score is outer to the penalized iteratively
re-weighted least squares scheme used to estimate the model given smoothing 
parameters.

These routines estimates a GAM given log smoothing paramaters, and evaluate
derivatives of the GCV and UBRE scores of the model with respect to the
log smoothing parameters. Calculation of exact derivatives is generally faster
than approximating them by finite differencing, as well as generally improving
the reliability of GCV/UBRE score minimization.

\code{gam.fit2} evaluates first derivatives, by accumulating them as P-IRLS
progresses. \code{gam.fit3} uses a more efficient approach in which the P-IRLS
is first run to convergence, and only then are the derivatives evaluated by a
separate iteration. \code{gam.fit3} can evaluate second as well as first derivatives.

Not normally called directly, but rather service routines for \code{\link{gam}}.
}
\usage{
gam.fit2(x, y, sp, S=list(),rS=list(),off, H=NULL, 
         weights = rep(1, nobs), start = NULL, etastart = NULL, 
         mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
         control = gam.control(), intercept = TRUE,deriv=TRUE,
         gamma=1,scale=1,pearson=FALSE,printWarn=TRUE,...)
gam.fit3(x, y, sp, S=list(),rS=list(),off, H=NULL, 
         weights = rep(1, nobs), start = NULL, etastart = NULL, 
         mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
         control = gam.control(), intercept = TRUE,deriv=2,use.svd=TRUE,
         gamma=1,scale=1,printWarn=TRUE,...) 
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x}{The model matrix for the GAM.}

\item{y}{The response variable.}

\item{sp}{The log smoothing parameters.}

\item{S}{A list of penalty matrices. Typically penalty matrices contain only a
smallish square sub-matrix which is non-zero: this is what is actually
stored. \code{off[i]} indicates which parameter is the first one penalized
by \code{S[[i]]}.}

\item{rS}{List of square roots of penalty matrices, each having as few columns
as possible, but as many rows as there are parameters.}

\item{off}{\code{off[i]} indicates which parameter \code{S[[i]][1,1]} relates to.}

\item{H}{The fixed penalty matrix for the model.}
 
\item{weights}{prior weights for fitting.}

\item{start}{optional starting parameter guesses.} 

\item{etastart}{optional starting values for the linear predictor.}

\item{mustart}{optional starting values for the mean.}

\item{offset}{the model offset}

\item{family}{the family - actually this routine would never be called with \code{gaussian()}} 

\item{control}{control list as returned from \code{\link{glm.control}}} 

\item{intercept}{does the model have and intercept, \code{TRUE} or
\code{FALSE}}

\item{deriv}{ Should derivatives of the GCV and UBRE scores be calculated?
            \code{TRUE} or \code{FALSE} for \code{gam.fit2}, but 0, 1 or 2,
            indicating the maximum order of differentiation to apply, for
            \code{gam.fit3}.}

\item{use.svd}{Should the algorithm use SVD (\code{TRUE}) or the cheaper QR
(\code{FALSE}) as the second matrix decomposition of the final
derivative/coefficient evaluation method? Only used by \code{gam.fit3}. }

\item{gamma}{The weight given to each degree of freedom in the GCV and UBRE
scores can be varied (usually increased) using this parameter.}

\item{scale}{The scale parameter - needed for the UBRE score.}

\item{pearson}{The GCV/UBRE score can be based either on the Pearson statistic
or the deviance. The latter is generally to be preferred, as it is less prone
to severe undersmoothing. Only used by \code{gam.fit2}.}

\item{printWarn}{Set to \code{FALSE} to suppress some warnings. Useful in
order to ensure that some warnings are only printed if they apply to the final
fitted model, rather than an intermediate used in optimization.}

\item{...}{Other arguments: ignored.}
}
\details{ These routines are basically \code{\link{glm.fit}} with some
modifications to allow (i) for quadratic penalties on the log likelihood;
(ii) derivatives of the model coefficients with respect to
log smoothing parameters to be obtained (by updating alongside the P-IRLS
iteration in the case of \code{gam.fit2}, and by an extra iteration in the
case of \code{gam.fit3}) and (iii) derivatives of the GAM GCV and UBRE scores to be
evaluated at convergence. 

In addition the routines apply step halving to any step that increases the
penalized deviance substantially.

The most costly parts of the calculations are performed by calls to compiled C
code (which in turn calls LAPACK routines) in place of the compiled code that
would usually perform least squares estimation on the working model in the
IRLS iteration. 

Estimation of smoothing parameters by optimizing GCV scores obtained at
convergence of the P-IRLS iteration was proposed by O'Sullivan et al. (1986),
and is here termed `outer' iteration. 

Note that use of non-standard families with this routine requires modification
of the families as described in \code{\link{fix.family.link}}.

}

\references{

O 'Sullivan, Yandall & Raynor (1986) Automatic smoothing of regression
functions in generalized linear models. J. Amer. Statist. Assoc. 81:96-103.

\url{http://www.maths.bath.ac.uk/~sw283/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}

The routine has been modified from \code{glm.fit} in R 2.0.1, written
by the R core (see \code{\link{glm.fit}} for further credits).
}


\seealso{\code{\link{gam.fit}},  \code{\link{gam}}, \code{\link{mgcv}}, \code{\link{magic}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






