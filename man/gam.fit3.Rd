\name{gam.fit3}
\alias{gam.fit3}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{P-IRLS GAM estimation with GCV \& UBRE/AIC derivative calculation}
\description{Estimation of GAM smoothing parameters is most stable if
optimization of the UBRE/AIC or GCV score is outer to the penalized iteratively
re-weighted least squares scheme used to estimate the model given smoothing 
parameters.

This routine estimates a GAM (any quadratically penalized GLM) given log 
smoothing paramaters, and evaluates derivatives of the GCV and UBRE/AIC scores 
of the model with respect to the
log smoothing parameters. Calculation of exact derivatives is generally faster
than approximating them by finite differencing, as well as generally improving
the reliability of GCV/UBRE/AIC score minimization.

The approach is to run the P-IRLS
to convergence, and only then to iterate for first and second derivatives. 

Not normally called directly, but rather service routines for \code{\link{gam}}.
}
\usage{

gam.fit3(x, y, sp, S=list(),rS=list(),off, H=NULL, 
         weights = rep(1, nobs), start = NULL, etastart = NULL, 
         mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
         control = gam.control(), intercept = TRUE,deriv=2,use.svd=TRUE,
         gamma=1,scale=1,printWarn=TRUE,...) 
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{x}{The model matrix for the GAM (or any penalized GLM).}

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

\item{deriv}{ Should derivatives of the GCV and UBRE/AIC scores be calculated?
             0, 1 or 2,
            indicating the maximum order of differentiation to apply.}

\item{use.svd}{Should the algorithm use SVD (\code{TRUE}) or the cheaper QR
(\code{FALSE}) as the second matrix decomposition of the final
derivative/coefficient evaluation method? }

\item{gamma}{The weight given to each degree of freedom in the GCV and UBRE
scores can be varied (usually increased) using this parameter.}

\item{scale}{The scale parameter - needed for the UBRE/AIC score.}

\item{printWarn}{Set to \code{FALSE} to suppress some warnings. Useful in
order to ensure that some warnings are only printed if they apply to the final
fitted model, rather than an intermediate used in optimization.}

\item{...}{Other arguments: ignored.}
}
\details{ This routine is basically \code{\link{glm.fit}} with some
modifications to allow (i) for quadratic penalties on the log likelihood;
(ii) derivatives of the model coefficients with respect to
log smoothing parameters to be obtained by an extra iteration  and 
(iii) derivatives of the GAM GCV and UBRE/AIC scores to be
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

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for
generalized additive models. J.R.Statist. Soc. B 70(3):495-518


\url{http://www.maths.bath.ac.uk/~sw283/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}

The routine has been modified from \code{glm.fit} in R 2.0.1, written
by the R core (see \code{\link{glm.fit}} for further credits).
}


\seealso{\code{\link{gam.fit}},  \code{\link{gam}}, \code{\link{mgcv}}, \code{\link{magic}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






