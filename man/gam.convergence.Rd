\name{gam.convergence}
\alias{gam.convergence}
\alias{gam.performance}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM convergence and performance issues}
\description{ When fitting GAMs there is a tradeoff between speed of 
fitting and probability of fit convergence. The fitting methods used 
by \code{\link{gam}} opt for certainty of convergence over speed of
fit. \code{\link{bam}} opts for speed.

 \code{\link{gam}} uses a nested iteration method (see \code{\link{gam.outer}}), in
 which each trial set of smoothing parameters proposed by an outer Newton algorithm
 require an inner Newton algorithm (penalized iteratively re-weighted least squares, PIRLS)
 to find the corresponding best fit model coefficients. Implicit differentiation is used to
 find the derivatives of the coefficients with respect to log smoothing parameters, so that the
 derivatives of the smoothness selection criterion can be obtained, as required by the outer iteration.
 This approach is less expensive than it at first appears, since excellent starting values for the inner
 iteration are available as soon as the smoothing parameters start to converge. See Wood (2011) and Wood, Pya and Saefken (2016).  
 

\code{\link{bam}} uses an alternative approach similar to `performance iteration' or `PQL'. A single PIRLS iteration is run to find the model coefficients. At each step this requires the estimation of a working penalized linear model. Smoothing parameter selection is applied directly to this working model at each step (as if it were a Gaussian additive model). This approach is more straightforward to code and in principle less costly than the nested approach. However it is not guaranteed to converge, since the smoothness selection criterion is changing at each iteration. It is sometimes possible for the algorithm to cycle around a small set of smoothing parameter, coefficient combinations without ever converging. \code{\link{bam}} includes some checks to limit this behaviour, and the further checks in the algorithm used by \code{bam(...,discrete=TRUE)} actually guarantee convergence in some cases, but in general guarantees are not possible. See Wood, Goude and Shaw (2015) and  Wood et al. (2017). 

\code{\link{gam}} when used with `general' families (such as \code{\link{multinom}} or \code{cox.ph}) can also use a potentially faster scheme based on the extended Fellner-Schall method (Wood and Fasiolo, 2017). This also operates with a single iteration and is not guaranteed to converge, theoretically. 

There are three things that you can try to speed up GAM fitting. (i) if you have large 
numbers of smoothing parameters in the generalized case, then try the \code{"bfgs"} method 
option in \code{\link{gam}} argument \code{optimizer}: this can be faster than the default. (ii) Try using
\code{\link{bam}}
(iii) For large datasets it may be worth changing
the smoothing basis to use \code{bs="cr"} (see \code{\link{s}} for details)
for 1-d smooths, and to use \code{\link{te}} smooths in place of
\code{\link{s}} smooths for smooths of more than one variable. This is because
the default thin plate regression spline basis \code{"tp"} is costly to set up
for large datasets.

If you have convergence problems, it's worth noting that a GAM is just a (penalized)
GLM and the IRLS scheme used to estimate GLMs is not guaranteed to
converge. Hence non convergence of a GAM may relate to a lack of stability in
the basic IRLS scheme. Therefore it is worth trying to establish whether the IRLS iterations
are capable of converging. To do this fit the problematic GAM with all smooth
terms specified with \code{fx=TRUE} so that the smoothing parameters are all
fixed at zero. If this `largest' model can converge then, then the maintainer 
would quite like to know about your problem! If it doesn't converge, then its
likely that your model is just too flexible for the IRLS process itself. Having tried
increasing \code{maxit} in \code{gam.control}, there are several other
possibilities for stabilizing the iteration. It is possible to try (i) setting lower bounds on the
smoothing parameters using the \code{min.sp} argument of \code{gam}: this may
or may not change the model being fitted; (ii)
reducing the flexibility of the model by reducing the basis dimensions
\code{k} in the specification of \code{s} and \code{te} model terms: this
obviously changes the model being fitted somewhat.

Usually, a major contributer to fitting difficulties is that the
model is a very poor description of the data.  
 
Please report convergence problems, especially if you there is no obvious pathology in the data/model that
suggests convergence should fail. 
}

%- maybe also `usage' for other objects documented here.
\references{

Key References on this implementation:

Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models (with discussion).
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

Wood, S.N., Goude, Y. & Shaw S. (2015) Generalized additive models for large datasets. Journal of the Royal Statistical Society, Series C 64(1): 139-155.

Wood, S.N., Li, Z., Shaddick, G. & Augustin N.H. (2017) Generalized additive models for gigadata: modelling the UK black smoke network daily data. Journal of the American Statistical Association.

Wood, S.N. and M. Fasiolo (2017) A generalized Fellner-Schall method for smoothing parameter optimization with application to Tweedie location, scale and shape models, Biometrics.

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press.
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
