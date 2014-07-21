\name{gam.scale}
\alias{gam.scale}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Scale parameter estimation in GAMs}
\description{Scale parameter estimation in \code{\link{gam}} depends on the type of \code{family}. For 
extended families then the RE/ML estimate is used. For conventional exponential families, estimated by the default 
outer iteration, the scale estimator can be controlled using argument \code{scale.est} in 
\code{\link{gam.control}}. The options are \code{"robust"}, \code{"pearson"} (default) or \code{"deviance"}. 
The Pearson estimator is the (weighted) sum of squares of the pearson residuals, divided by the 
effective residual degrees of freedom. The deviance estimator simply substitutes deviance residuals for Pearson residuals. 

Usually the Pearson estimator is recommended for GLMs, since it is asymptotically unbiased. However, it can also be unstable at finite sample sizes, if a few Pearson residuals are very large. For example, a very low Poisson mean with a non zero count can give a huge Pearson residual, even though the deviance residual is much more modest. The (\code{"robust"}) estimator in \code{\link{gam}} is an attempt to mitigate against the bias of the deviance estimator and the instability of the Pearson estimator. The estimator simply scales the deviance residuals by the median ratio of pearson to deviance residuals, separately for positive and negative residuals (given the difference in skew of the two types). The estimator is forced to lie between the Pearson and deviance estimators. All three options coincide for Gaussian data. A warning is generated if the Pearson scale estimate is more than twice the robust version. In this case you consider changing to the robust estimator, especially if the Pearson residuals contain extreme outliers.

For performance iteration the Pearson estimator is always used.

\code{\link{gamm}} uses the estimate of the scale parameter from the underlying call to \code{lme}. \code{\link{bam}} uses the REML estimator if the method is \code{"fREML"}. Otherwise the estimator is a Pearson estimator.

}


\author{ Simon N. Wood \email{simon.wood@r-project.org} } 

\seealso{  \code{\link{gam.control} } }

\keyword{models} \keyword{smooth} \keyword{regression} %-- one or more ..






