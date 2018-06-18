\name{k.check}
\alias{k.check}
\title{Checking smooth basis dimension }
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and runs 
diagnostic tests of whether the basis dimension choises are adequate. 
}
\usage{
k.check(b, subsample=5000, n.rep=400)
}
\arguments{
  \item{b}{a fitted \code{gam} object as produced by \code{\link{gam}()}.}
  \item{subsample}{above this number of data, testing uses a random sub-sample of data of this size.}
  \item{n.rep}{how many re-shuffles to do to get p-value for k testing.}
}

\value{A matrix contaning the output of the tests described above.}

\details{ 
The test of whether the basis dimension for a smooth is adequate (Wood, 2017, section 5.9) is based on computing an estimate of the residual variance based on differencing residuals that are near neighbours according to the (numeric) covariates of the smooth. This estimate divided by the residual variance is the \code{k-index} reported. The further below 1 this is, the more likely it is that there is missed pattern left in the residuals. The \code{p-value} is computed by simulation: the residuals are randomly re-shuffled \code{n.rep} times to obtain the null distribution of the differencing variance estimator, if there is no pattern in the residuals. For models fitted to more than \code{subsample} data, the tests are based of \code{subsample} randomly sampled data. Low p-values may indicate that the basis dimension, \code{k}, has been set too low, especially if the reported \code{edf} is close to \code{k\'}, the maximum possible EDF for the term. Note the disconcerting fact that if the test statistic itself is based on random resampling and the null is true, then the associated p-values will of course vary widely from one replicate to the next. Currently smooths of factor variables are not supported and will give an \code{NA} p-value.

Doubling a suspect \code{k} and re-fitting is sensible: if the reported \code{edf} increases substantially then you may have been missing something in the first fit. Of course p-values can be low for reasons other than a too low \code{k}. See \code{\link{choose.k}} for fuller discussion.
 
}

\references{

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press.

\url{http://www.maths.bris.ac.uk/~sw15190/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\seealso{ \code{\link{choose.k}},  \code{\link{gam}}, \code{\link{gam.check}}}

\examples{
library(mgcv)
set.seed(0)
dat <- gamSim(1,n=200)
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
plot(b,pages=1)
k.check(b)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...









