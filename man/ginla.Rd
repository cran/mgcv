\name{ginla}
\alias{ginla}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM Integrated Nested Laplace Approximation Newton Enhanced}
\description{Apply Integrated Nested Laplace Approximation (INLA, Rue et al. 2009) to models estimable by \code{\link{gam}} or \code{\link{bam}}, using the INLA variant described in Wood (2019). Produces marginal posterior densities for each coefficient, selected coefficients or linear transformations of the coefficient vector. 
}
\usage{
ginla(G,A=NULL,nk=16,nb=100,J=1,interactive=FALSE,int=0,approx=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{G}{A pre-fit gam object, as produced by \code{gam(...,fit=FALSE)} or \code{bam(...,discrete=TRUE,fit=FALSE)}.}
\item{A}{Either a matrix of transforms of the coefficients that are of interest, or an array of indices of the parameters of interest. If \code{NULL} then distributions are produced for all coefficients.}
\item{nk}{Number of values of each coefficient at which to evaluate its log marginal posterior density. These points are then spline interpolated.}
\item{nb}{Number of points at which to evaluate posterior density of coefficients for returning as a gridded function.}
\item{J}{How many determinant updating steps to take in the log determinant approximation step. Not recommended to increase this. }
\item{interactive}{If this is \code{>0} or \code{TRUE} then every approximate posterior is plotted in red, overlaid on the simple Gaussian approximate posterior. If \code{2} then waits for user to press return between each plot. Useful for judging whether anything is gained by using INLA approach. }
\item{int}{0 to skip integration and just use the posterior modal smoothing parameter. >0 for integration using the CCD approach proposed in Rue et al. (2009).}
\item{approx}{0 for full approximation; 1 to update Hessian, but use approximate modes; 2 as 1 and assume constant Hessian. See details.}
}

\value{
A list with elements \code{beta} and \code{density}, both of which are matrices. Each row relates to one coefficient (or linear coefficient combination) of interest. Both matrices have \code{nb} columns. If \code{int!=0} then a further element \code{reml} gives the integration weights used in the CCD integration, with the central point weight given first.  
}

\details{Let \eqn{\beta}{b}, \eqn{\theta}{h} and \eqn{y}{y} denote the model coefficients, hyperparameters/smoothing parameters and response data, respectively. In principle, INLA employs Laplace approximations for \eqn{\pi(\beta_i|\theta,y)}{p(b_i|h,y)} and \eqn{\pi(\theta|y)}{p(h|y)} and then obtains the marginal posterior distribution \eqn{\pi(\beta_i|y)}{p(b_i|y)} by intergrating the approximations to \eqn{\pi(\beta_i|\theta,y)\pi(\theta|y)}{p(b_i|h,y)p(h|y)} w.r.t \eqn{\theta}{h} (marginals for the hyperparameters can also be produced). In practice the Laplace approximation for \eqn{\pi(\beta_i|\theta,y)}{p(b_i|h,y)} is too expensive to compute for each \eqn{\beta_i}{b_i} and must itself be approximated. To this end, there are two quantities that have to be computed: the posterior mode \eqn{\beta^*|\beta_i}{b*|b_i} and the determinant of the Hessian of the joint log density \eqn{\log \pi(\beta,\theta,y)}{log p(b,h,y)} w.r.t. \eqn{\beta}{b} at the mode. Rue et al. (2009) originally approximated the posterior conditional mode by the conditional mode implied by a simple Gaussian approximation to the posterior \eqn{\pi(\beta|y)}{p(b|y)}. They then approximated the log determinant of the Hessian as a function of \eqn{\beta_i}{b_i} using a first order Taylor expansion, which is cheap to compute for the sparse model representaiton that they use, but not when using the dense low rank basis expansions used by \code{\link{gam}}. They also offer a more expensive alternative approximation based on computing the log determiannt with respect only to those elements of \eqn{\beta}{b} with sufficiently high correlation with \eqn{\beta_i}{b_i} according to the simple Gaussian posterior approximation: efficiency again seems to rest on sparsity. Wood (2018) suggests computing the required posterior modes exactly, and basing the log determinant approximation on a BFGS update of the Hessian at the unconditional model. The latter is efficient with or without sparsity, whereas the former is a `for free' improvement. Both steps are efficient because it is cheap to obtain the Cholesky factor of \eqn{H[-i,-i]}{H[-i,-i]} from that of \eqn{H}{H} - see \code{\link{choldrop}}. This is the approach taken by this routine.

The \code{approx} argument allows two further approximations to speed up computations. For \code{approx==1} the exact posterior conditional modes are not used, but instead the conditional modes implied by the simple Gaussian posterior approximation. For \code{approx==2} the same approximation is used for the modes and the Hessian is assumed constant. The latter is quite fast as no log joint density gradient evaluations are required.  


Note that for many models the INLA estimates are very close to the usual Gaussian approximation to the posterior, the \code{interactive} argument is useful for investigating this issue.

\code{\link{bam}} models are only supported with the \code{disrete=TRUE} option. The \code{discrete=FALSE} approach would be too inefficient. AR1 models are not supported (related arguments are simply ignored). 
}

\references{
Rue, H, Martino, S. & Chopin, N. (2009) Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations (with discussion). Journal of the Royal Statistical Society, Series B.  71: 319-392.

Wood (2019) Simplified Integrated Laplace Approximation. In press Biometrika. 
}

\section{WARNINGS}{
This routine is still somewhat experimental, so details are liable to change. Also currently not all steps are optimally efficient.

The routine is written for relatively expert users.

\code{ginla} is not designed to deal with rank deficient models.

}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
  require(mgcv); require(MASS)

  ## example using a scale location model for the motorcycle data. A simple plotting
  ## routine is produced first...

  plot.inla <- function(x,inla,k=1,levels=c(.025,.1,.5,.9,.975),
               lcol = c(2,4,4,4,2),lwd = c(1,1,2,1,1),lty=c(1,1,1,1,1),
	       xlab="x",ylab="y",cex.lab=1.5) {
    ## a simple effect plotter, when distributions of function values of
    ## 1D smooths have been computed
    require(splines)
    p <- length(x) 
    betaq <- matrix(0,length(levels),p) ## storage for beta quantiles 
    for (i in 1:p) { ## work through x and betas
      j <- i + k - 1
      p <- cumsum(inla$density[j,])*(inla$beta[j,2]-inla$beta[j,1])
      ## getting quantiles of function values...
      betaq[,i] <- approx(p,y=inla$beta[j,],levels)$y
    }
    xg <- seq(min(x),max(x),length=200)
    ylim <- range(betaq)
    ylim <- 1.1*(ylim-mean(ylim))+mean(ylim)
    for (j in 1:length(levels)) { ## plot the quantiles
      din <- interpSpline(x,betaq[j,])
      if (j==1) {
        plot(xg,predict(din,xg)$y,ylim=ylim,type="l",col=lcol[j],
             xlab=xlab,ylab=ylab,lwd=lwd[j],cex.lab=1.5,lty=lty[j])
      } else lines(xg,predict(din,xg)$y,col=lcol[j],lwd=lwd[j],lty=lty[j])
    }
  } ## plot.inla

  ## set up the model with a `gam' call...

  G <- gam(list(accel~s(times,k=20,bs="ad"),~s(times)),
            data=mcycle,family=gaulss(),fit=FALSE)
  b <- gam(G=G,method="REML") ## regular GAM fit for comparison

  ## Now use ginla to get posteriors of estimated effect values
  ## at evenly spaced times. Create A matrix for this...
  
  rat <- range(mcycle$times)
  pd0 <- data.frame(times=seq(rat[1],rat[2],length=20))
  X0 <- predict(b,newdata=pd0,type="lpmatrix")
  X0[,21:30] <- 0
  pd1 <- data.frame(times=seq(rat[1],rat[2],length=10))
  X1 <- predict(b,newdata=pd1,type="lpmatrix")
  X1[,1:20] <- 0
  A <- rbind(X0,X1) ## A maps coefs to required function values

  ## call ginla. Set int to 1 for integrated version.
  ## Set interactive = 1 or 2 to plot marginal posterior distributions
  ## (red) and simple Gaussian approximation (black).
 
  inla <- ginla(G,A,int=0)

  par(mfrow=c(1,2),mar=c(5,5,1,1))
  fv <- predict(b,se=TRUE) ## usual Gaussian approximation, for comparison

  ## plot inla mean smooth effect...
  plot.inla(pd0$times,inla,k=1,xlab="time",ylab=expression(f[1](time))) 

  ## overlay simple Gaussian equivalent (in grey) ...
  points(mcycle$times,mcycle$accel,col="grey")
  lines(mcycle$times,fv$fit[,1],col="grey",lwd=2)
  lines(mcycle$times,fv$fit[,1]+2*fv$se.fit[,1],lty=2,col="grey",lwd=2)
  lines(mcycle$times,fv$fit[,1]-2*fv$se.fit[,1],lty=2,col="grey",lwd=2)

  ## same for log sd smooth...
  plot.inla(pd1$times,inla,k=21,xlab="time",ylab=expression(f[2](time)))
  lines(mcycle$times,fv$fit[,2],col="grey",lwd=2)
  lines(mcycle$times,fv$fit[,2]+2*fv$se.fit[,2],col="grey",lty=2,lwd=2)
  lines(mcycle$times,fv$fit[,2]-2*fv$se.fit[,2],col="grey",lty=2,lwd=2)

  ## ... notice some real differences for the log sd smooth, especially
  ## at the lower and upper ends of the time interval.
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









