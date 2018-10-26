\name{smooth2random}
\alias{smooth2random}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert a smooth to a form suitable for estimating as random effect}
\description{A generic function for converting \code{mgcv} smooth objects to forms suitable for estimation as random effects by e.g. \code{lme}. Exported mostly for use by other package developers. 
}
\usage{
smooth2random(object,vnames,type=1)

%\method{summary}{gam}(object, dispersion=NULL, freq=FALSE, p.type = 0, ...)

}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{an \code{mgcv} smooth object.}

\item{vnames}{a vector of names to avoid as dummy variable names in the random effects form.} 

\item{type}{\code{1} for \code{lme}, otherwise \code{lmer}.}

}

\details{There is a duality between smooths and random effects which means that smooths can be estimated using mixed modelling software. This function converts standard \code{mgcv} smooth objects to forms suitable for estimation by \code{lme}, for example. A service routine for \code{\link{gamm}} exported for use by package developers. See examples for creating prediction matrices for new data, corresponding to the random and fixed effect matrices returned when \code{type=2}. 
}

\value{A list.

\item{rand}{ a list of random effects, including grouping factors, and 
       a fixed effects matrix. Grouping factors, model matrix and model
       matrix name attached as attributes, to each element. Alternatively, for \code{type=2}
       a list of random effect model matrices, each corresponding to an i.i.d. Gaussian
       random effect with a single variance component.}

\item{trans.D}{A vector, trans.D, that transforms coefs, in order [rand1, rand2,... fix] back to original parameterization. If null, then taken as vector of ones. \code{b.original = trans.U \%*\% (trans.D*b.fit)}.}

\item{trans.U}{A matrix, trans.U, that transforms coefs, in order [rand1, rand2,... fix] back to original parameterization. If null, then not needed. If null then taken as identity.}

\item{Xf}{A matrix for the fixed effects, if any.}

\item{fixed}{\code{TRUE/FALSE}, indicating if term was unpenalized or not. If unpenalized then other stuff may not be returned (it's not a random effect).}

\item{rind}{an index vector such that if br is the vector of 
random coefficients for the term, br[rind] is the coefs in 
order for this term. }

\item{pen.ind}{index of which penalty penalizes each coefficient: 0 for unpenalized.}
}


\references{

Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman
and Hall/CRC Press.

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}.}


\seealso{  \code{\link{gamm}} }

\examples{
## Simple type 1 'lme' style...
library(mgcv)
x <- runif(30)
sm <- smoothCon(s(x),data.frame(x=x))[[1]]
smooth2random(sm,"")

## Now type 2 'lme4' style...
z <- runif(30)
dat <- data.frame(x=x,z=z)
sm <- smoothCon(t2(x,z),dat)[[1]]
re <- smooth2random(sm,"",2)
str(re)

## For prediction after fitting we might transform parameters back to
## original parameterization using 'rind', 'trans.D' and 'trans.U',
## and call PredictMat(sm,newdata) to get the prediction matrix to
## multiply these transformed parameters by.
## Alternatively we could obtain fixed and random effect Prediction
## matrices corresponding to the results from smooth2random, which
## can be used with the fit parameters without transforming them.
## The following shows how...

s2rPred <- function(sm,re,data) {
## Function to aid prediction from smooths represented as type==2
## random effects. re must be the result of smooth2random(sm,...,type=2).
  X <- PredictMat(sm,data)   ## get prediction matrix for new data
  ## transform to r.e. parameterization
  if (!is.null(re$trans.U)) X <- X\%*\%re$trans.U
  X <- t(t(X)*re$trans.D)
  ## re-order columns according to random effect re-ordering...
  X[,re$rind] <- X[,re$pen.ind!=0] 
  ## re-order penalization index in same way  
  pen.ind <- re$pen.ind; pen.ind[re$rind] <- pen.ind[pen.ind>0]
  ## start return object...
  r <- list(rand=list(),Xf=X[,which(re$pen.ind==0),drop=FALSE])
  for (i in 1:length(re$rand)) { ## loop over random effect matrices
    r$rand[[i]] <- X[,which(pen.ind==i),drop=FALSE]
    attr(r$rand[[i]],"s.label") <- attr(re$rand[[i]],"s.label")
  }
  names(r$rand) <- names(re$rand)
  r
} ## s2rPred

## use function to obtain prediction random and fixed effect matrices
## for first 10 elements of 'dat'. Then confirm that these match the
## first 10 rows of the original model matrices, as they should...

r <- s2rPred(sm,re,dat[1:10,])
range(r$Xf-re$Xf[1:10,])
range(r$rand[[1]]-re$rand[[1]][1:10,])

}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






