\name{mgcv.parallel}
\alias{mgcv.parallel}
\title{Parallel computation in mgcv.
}
\description{
\code{mgcv} can make some use of multiple cores or a cluster. Function \code{\link{bam}} uses the 
facilities provided in the \link[parallel]{parallel} package for this purpose. See examples in 
\code{\link{bam}}. 

Function \code{\link{gam}} can use parallel threads on a (shared memory) multi-core 
machine via \code{openMP} (where this is supported). To do this, set the desired number of threads by setting \code{nthreads} to the number of cores to use, in the \code{control} argument of \code{\link{gam}}. Note that, for the most part, only the dominant \eqn{O(np^2)}{O(np^2)} steps are parallelized (n is number of data, p number of parameters). For additive Gaussian models estimated by GCV, the speed up can be disappointing as these employ an \eqn{O(p^3)}{O(p^3)} SVD step that can also have substantial cost in practice. The \eqn{O(np^2)}{O(np^3)} QR decomposition steps are only parallelized if this would be worthwhile on a machine in which each parallel thread runs as quickly as a single thread. 

\code{\link{magic}} can also use multiple cores, but the same comments apply as for the GCV Gaussian additive model. 

If \code{control$nthreads} is set to more than the number of cores detected, then only the number of detected cores is used. Note that using virtual cores usually gives very little speed up, and can even slow computations slightly. For example, many Intel processors reporting 4 cores actually have 2 physical cores, each with 2 virtual cores, so using 2 threads gives a marked increase in speed, while using 4 threads makes little extra difference. 

Because the computational burden in \code{mgcv} is all in the linear algebra, then parallel computation may provide reduced or no benefit with a tuned BLAS. This is particularly the case if you are using a multi threaded BLAS, but a BLAS that is tuned to make efficient use of a particular cache size may also experience loss of performance if threads have to share the cache.  

}

\author{
Simon Wood <simon.wood@r-project.org>
}
\references{
\url{https://computing.llnl.gov/tutorials/openMP/}
}

\keyword{ package }
\keyword{models} 
\keyword{smooth} 
\keyword{regression}


\examples{
## illustration of multi-threading with gam...

require(mgcv);set.seed(9)
dat <- gamSim(1,n=2000,dist="poisson",scale=.1)
k <- 12;bs <- "cr";ctrl <- list(nthreads=2)

system.time(b1<-gam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+
            s(x3,bs=bs,k=k),family=poisson,data=dat,method="REML"))[3]

system.time(b2<-gam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+
s(x3,bs=bs,k=k),family=poisson,data=dat,method="REML",control=ctrl))[3]

## Poisson example on a cluster with 'bam'. 
## Note that there is some overhead in initializing the 
## computation on the cluster, associated with loading 
## the Matrix package on each node. For this reason the 
## sample sizes here are very small to keep CRAN happy, but at
## this low sample size you see little advantage of parallel computation.

k <- 15
dat <- gamSim(1,n=8000,dist="poisson",scale=.1)
require(parallel)  
nc <- 2   ## cluster size, set for example portability
if (detectCores()>1) { ## no point otherwise
  cl <- makeCluster(nc) 
  ## could also use makeForkCluster, but read warnings first!
} else cl <- NULL
  
system.time(b3 <- bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)
            ,data=dat,family=poisson(),chunk.size=5000,cluster=cl))

fv <- predict(b3,cluster=cl) ## parallel prediction

if (!is.null(cl)) stopCluster(cl)
b3

}
