\name{family.mgcv}
\alias{family.mgcv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Distribution families in mgcv}
\description{As well as the standard families documented in \code{\link{family}} (see also \code{\link{glm}}) which can be used with functions \code{\link{gam}}, \code{\link{bam}} and \code{\link{gamm}}, \code{mgcv} also supplies some extra families, most of which are currently only usable with \code{\link{gam}}, although some can also be used with \code{\link{bam}}. These are described here.
}


\details{
The following families are in the exponential family given the value of a single parameter. They are usable with all modelling functions.
\itemize{
\item \code{\link{Tweedie}} An exponential family distribution for which the variance of the response is given by the mean response to the power \code{p}. 
\code{p} is in (1,2) and must be supplied. Alternatively, see \code{\link{tw}} to estimate \code{p} (\code{gam} only).
\item \code{\link{negbin}} The negative binomial. Alternatively see \code{\link{nb}} to estimate the \code{theta} parameter of the negative binomial (\code{gam} only). 

}


The following families are for regression type models dependent on a single linear predictor, and with a log likelihood
which is a sum of independent terms, each coprresponding to a single response observation. Usable with \code{\link{gam}}, with smoothing parameter estimation by \code{"REML"} or \code{"ML"} (the latter does not integrate the unpenalized and parameteric effects out of the marginal likelihood optimized for the smoothing parameters). Also usable with \code{\link{bam}}.
\itemize{
\item \code{\link{ocat}} for ordered categorical data.
\item \code{\link{tw}} for Tweedie distributed data, when the power parameter relating the variance to the mean is to be estimated.
\item \code{\link{nb}} for negative binomial data when the \code{theta} parameter is to be estimated.
\item \code{\link{betar}} for proportions data on (0,1) when the binomial is not appropriate.
\item \code{\link{scat}} scaled t for heavy tailed data that would otherwise be modelled as Gaussian.
\item \code{\link{ziP}} for zero inflated Poisson data, when the zero inflation rate depends simply on the Poisson mean.
} %% end itemize 

The following families implement more general model classes. Usable only with \code{\link{gam}} and only with REML smoothing parameter estimation.
\itemize{
\item \code{\link{cox.ph}} the Cox Proportional Hazards model for survival data.
\item \code{\link{gaulss}} a Gaussian location-scale model where the mean and the standard deviation are both modelled using smooth linear predictors.
\item \code{\link{gevlss}} a generalized extreme value (GEV) model where the location, scale and shape parameters are each modelled using a linear predictor. 
\item \code{\link{ziplss}} a `two-stage' zero inflated Poisson model, in which 'potential-presence' is modelled with one linear predictor, and Poisson mean abundance
given potential presence is modelled with a second linear predictor. 
\item \code{\link{mvn}}: multivariate normal additive models.  
\item \code{\link{multinom}}: multinomial logistic regression, for unordered categorical responses.  
 } %% end itemize 
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood (s.wood@r-project.org) & Natalya Pya
}
\references{Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\url{http://dx.doi.org/10.1080/01621459.2016.1180986}
}


\keyword{models} \keyword{regression}%-- one or more ..


