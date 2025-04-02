\name{cpois}
\alias{cpois}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM censored Poisson family}
\description{Family for use with \code{\link{gam}} or \code{\link{bam}}, implementing regression for censored
Poisson data. Observations may be left, interval or right censored or uncensored.
}

\usage{
cpois(link="log")
}
\arguments{

\item{link}{The link function: \code{"identity"}, \code{"log"} or \code{"sqrt"}.}

}
\value{
An object of class \code{extended.family}.
}

\details{If the family is used with a vector response, then it is assumed that there is no censoring, and a regular Poisson regression results. If there is censoring then the response should be supplied as a two column matrix. The first column is always numeric. Entries in the second column are as follows.
\itemize{
\item If an entry is identical to the corresponding first column entry, then it is an uncensored observation.
\item If an entry is numeric and different to the first column entry then there is interval censoring. The first column entry is the lower interval limit and the second column entry is the upper interval limit (both should be non-integer). \eqn{y}{y} is only known to be between these limits.
\item If the second column entry is \code{-Inf} then the observation is left censored at the value of the entry (non-integer) in the first column. It is only known that \eqn{y}{y} is below the first column value.
\item If the second column entry is \code{Inf} then the observation is right censored at the value (non-integer) of the entry in the first column. It is only known that \eqn{y}{y} is above the first column value.
}
Any mixture of censored and uncensored data is allowed, but be aware that data consisting only of right and/or left censored data contain very little information. It is strongly recommended to use non-integer values for censoring limits, to avoid any possibility of ambiguity. For example if \eqn{y}{y} is known to be 3 or above, set the lower censoring limit to 2.5, or any other value in \eqn{(2,3)}{(2,3)}.
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\doi{10.1080/01621459.2016.1180986}
}


\examples{
library(mgcv)
set.seed(6); n <- 2000
dat <- gamSim(1,n=n,dist="poisson",scale=.1) ## simulate Poi data

## illustration that cpois an poisson give same results if there
## is no censoring...

b0 <- gam(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
             s(x3,bs="cr"),family=poisson,data=dat,method="REML")
plot(b0,pages=1,scheme=2)

b1 <- gam(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
            s(x3,bs="cr"),family=cpois,data=dat) 
plot(b1,pages=1,scheme=2)

b0;b1

## Now censor some observations...
dat1 <- dat
m <- 300
y <- matrix(dat$y,n,ncol=2) ## response matrix
ii <- sample(n,3*m) ## censor these
for (i in 1:m) { ## right, left, interval...
  j <- ii[i]; if (y[j,1] > 4.5) y[j,] <- c(4.5,Inf)
  j <- ii[m+i]; if (y[j,1] < 2.5) y[j,] <- c(2.5,-Inf)
  j <- ii[2*m+i];if (y[j,1] > 1.5 & y[j,1]< 5.5) y[j,] <- c(1.5,5.5)
}
n - sum(y[,1]==y[,2]) ## number of censored obs
dat1$y <- y

## now fit model with censoring...
b2 <- gam(y~s(x0,bs="cr")+s(x1,bs="cr")+s(x2,bs="cr")+
          s(x3,bs="cr"),family=cpois,data=dat1) 
plot(b2,pages=1,scheme=2);b2

}
\keyword{models} \keyword{regression}%-- one or more ..


