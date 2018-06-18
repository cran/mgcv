\name{logLik.gam}
\alias{logLik.gam}
\alias{AIC.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{AIC and Log likelihood for a fitted GAM}
\description{ Function to extract the log-likelihood for a fitted \code{gam}
model (note that the models are usually fitted by penalized likelihood maximization). 
Used by \code{\link{AIC}}. See details for more information on AIC computation.  
}
\usage{
\method{logLik}{gam}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{ fitted model objects of class \code{gam} as produced by \code{gam()}.}
\item{...}{un-used in this case} 

}
\details{ Modification of \code{logLik.glm} which corrects the degrees of
freedom for use with \code{gam} objects. 

The function is provided so that \code{\link{AIC}} functions correctly with
\code{gam} objects, and uses the appropriate degrees of freedom (accounting
for penalization). See e.g. Wood, Pya and Saefken (2016) for a derivation of
an appropriate AIC.

There are two possibile AIC's that might be considered for use with GAMs. Marginal
AIC is based on the marginal likelihood of the GAM, that is the likelihood based on
treating penalized (e.g. spline) coefficients as random and integrating them out. The
degrees of freedom is then the number of smoothing/variance parameters + the number
of fixed effects. The problem with Marginal AIC is that marginal likelihood
underestimates variance components/oversmooths, so that the approach favours simpler models
excessively (substituting REML does not work, because REML is not comparable between models
with different unpenalized/fixed components). Conditional AIC uses the likelihood of all
the model coefficients, evaluated at the penalized MLE. The degrees of freedom to use then
is the effective degrees of freedom for the model. However, Greven and Kneib (2010) show
that the neglect of smoothing parameter uncertainty can lead to this conditional AIC being
excessively likely to select larger models. Wood, Pya and Saefken (2016) propose a simple
correction to the effective degrees of freedom to fix this problem. \code{mgcv} applies this
correction whenever possible: that is when using \code{ML} or \code{REML} smoothing parameter
selection with \code{\link{gam}} or \code{\link{bam}}. The correction
is not computable when using the Extended Fellner Schall or BFGS optimizer (since the correction requires
an estimate of the covariance matrix of the log smoothing parameters). 


}

\value{ Standard \code{logLik} object: see \code{\link{logLik}}.
}


\references{

Greven, S., and Kneib, T. (2010), On the Behaviour of Marginal and
Conditional AIC in Linear Mixed Models, Biometrika, 97, 773-789.

Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models (with discussion).
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}

Wood S.N. (2017) Generalized Additive Models: An Introduction with R
(2nd edition). Chapman and Hall/CRC Press.


}

\author{ Simon N. Wood \email{simon.wood@r-project.org} based directly on \code{logLik.glm}}


\seealso{  \code{\link{AIC}}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






