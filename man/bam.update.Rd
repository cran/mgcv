\name{bam.update}
\alias{bam.update}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Update a strictly additive bam model for new data.}

\description{ Strictly additive models fitted by \code{\link{bam}} can be efficiently updated as new data becomes available,
by simply updateing the QR decomposition on which estimation is based, and re-optimizing the smoothing parameters, starting
from the previous estimates. This routine implements this.
}
\usage{
bam.update(b,data,chunk.size=10000)
}
%- maybe also `usage' for other objects documented here.

\arguments{ 
\item{b}{ A \code{gam} object fitted by \code{\link{bam}} and representing a strictly additive model 
(i.e. \code{gaussian} errors, \code{identity} link).}
\item{data}{Extra data to augment the original data used to obtain \code{b}. Must include a column of weights if the 
            original fit was weighted.}
\item{chunk.size}{size of subsets of data to process in one go when getting fitted values.}
}

\value{ 
An object of class \code{"gam"} as described in \code{\link{gamObject}}.
}


\details{ \code{bam.update} updates the QR decomposition of the (weighted) model matrix of the GAM represented by \code{b} to take 
account of the new data. The orthogonal factor multiplied by the response vector is also updated. Given these updates the model 
and smoothing parameters can be re-estimated, as if the whole dataset (original and the new data) had been fitted in one go. 

Note that there may be small numerical differences in fit between fitting the data all at once, and fitting in stages by updating, if the 
smoothing bases used have any of their detials set with reference to the data (e.g. default knot locations).

}


\references{


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}


}


\seealso{\code{\link{mgcv-package}}, \code{\link{bam}}
 }

\examples{
library(mgcv)
## following is not *very* large, for obvious reasons...
set.seed(8)

dat <- gamSim(1,n=5000,dist="normal",scale=5)
dat[c(50,13,3000,3005,3100),]<- NA
dat1 <- dat[4001:5000,]
dat0 <- dat[1:4000,]
bs <- "ps";k <- 20
method <- "GCV.Cp"
b <- bam(y ~ s(x0,bs=bs,k=k)+s(x1,bs=bs,k=k)+s(x2,bs=bs,k=k)+
           s(x3,bs=bs,k=k),data=dat0,method=method)

b1 <- bam.update(b,dat1)

b2 <- bam.update(bam.update(b,dat1[1:500,]),dat1[501:1000,])
 
b3 <- bam(y ~ s(x0,bs=bs,k=k)+s(x1,bs=bs,k=k)+s(x2,bs=bs,k=k)+
           s(x3,bs=bs,k=k),data=dat,method=method)
b1;b2;b3

summary(b1);summary(b2);summary(b3)

}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
\concept{Varying coefficient model}
\concept{Functional linear model}
\concept{Penalized GLM}
\concept{Generalized Additive Model}
\concept{Penalized regression}
\concept{Spline smoothing}
\concept{Penalized regression spline}
\concept{Generalized Cross Validation}
\concept{Smoothing parameter selection}
\concept{tensor product smoothing}
\concept{thin plate spline}
\concept{P-spline}
\concept{Generalized ridge regression}






