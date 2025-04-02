\name{clog}
\alias{clog}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM censored logistic distribution family for log-logistic AFT models}
\description{Family for use with \code{\link{gam}} or \code{\link{bam}}, implementing regression for censored
logistic data. If \eqn{y}{y} is the response with mean \eqn{\mu}{m} and scale parameter \eqn{s=w^{-1/2}\exp(\theta)}{s=w^{-1/2}exp(theta)},
then \eqn{y}{y} has p.d.f.
\deqn{\frac{\exp\{-(y-\mu)/s\}}{s[1+\exp\{-(y-\mu)/s\}]^2}.}{exp\{-(y-m)/s\}/(s[1+ exp\{-(y-m)/s\}]^2).}
 \eqn{\theta}{theta} is a single scalar for all observations. Observations may be left, interval or right censored or uncensored.

Useful for log-logistic accelerated failure time (AFT) models, for example. 
}

\usage{
clog(theta=NULL,link="identity")
}
\arguments{
\item{theta}{ log scale parameter. If supplied and positive then taken as a fixed value of standard deviation (not its log). If supplied and negative taken as negative of initial value for standard deviation (not its log).}

\item{link}{The link function: \code{"identity"}, \code{"log"} or \code{"sqrt"}.}

}
\value{
An object of class \code{extended.family}.
}

\details{If the family is used with a vector response, then it is assumed that there is no censoring, and a regular Gaussian regression results. If there is censoring then the response should be supplied as a two column matrix. The first column is always numeric. Entries in the second column are as follows.
\itemize{
\item If an entry is identical to the corresponding first column entry, then it is an uncensored observation.
\item If an entry is numeric and different to the first column entry then there is interval censoring. The first column entry is the lower interval limit and the second column entry is the upper interval limit. \eqn{y}{y} is only known to be between these limits.
\item If the second column entry is \code{-Inf} then the observation is left censored at the value of the entry in the first column. It is only known that \eqn{y}{y} is less than or equal to the first column value.
\item If the second column entry is \code{Inf} then the observation is right censored at the value of the entry in the first column. It is only known that \eqn{y}{y} is greater than or equal to the first column value.
}
Any mixture of censored and uncensored data is allowed, but be aware that data consisting only of right and/or left censored data contain very little information.
}

%- maybe also `usage' for other objects documented here.

\author{ Chris Shen
}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\doi{10.1080/01621459.2016.1180986}
}


\examples{
library(mgcv)

#######################################################
## AFT model example for colon cancer survivial data...
#######################################################

library(survival) ## for data
col1 <- colon[colon$etype==1,] ## concentrate on single event
col1$differ <- as.factor(col1$differ)
col1$sex <- as.factor(col1$sex)

## set up the AFT response... 
logt <- cbind(log(col1$time),log(col1$time))
logt[col1$status==0,2] <- Inf ## right censoring
col1$logt <- -logt ## -ve conventional for AFT versus Cox PH comparison

## fit the model...
b <- gam(logt~s(age,by=sex)+sex+s(nodes)+perfor+rx+obstruct+adhere,
         family=clog(),data=col1)
plot(b,pages=1)	 
## ... compare this to ?cox.ph
}
\keyword{models} \keyword{regression}%-- one or more ..


