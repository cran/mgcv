\name{NCV}
\alias{NCV}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Neighbourhood Cross Validation}
\description{
NCV estimates smoothing parameters by optimizing the average ability of a model to predict subsets of data when subsets of data are omitted from fitting. Usually the predicted subset is a subset of the omitted subset. If both subsets are the same single datapoint, and the average is over all datapoints, then NCV is leave-one-out cross validation. QNCV is a quadratic approximation to NCV, guaranteed finite for any family link combination.

In detail, suppose that a model is estimated by minimizing a penalized loss
\deqn{\sum_i D(y_i,\theta_i) + \sum_j \lambda_j \beta^{\sf T} {S}_j \beta }{sum_i D(y_i,t_i) + sum_j v_j b'S_jb}
where \eqn{D}{D} is a loss (such as a negative log likelihood), dependent on response \eqn{y_i}{y_i} and parameter vector \eqn{\theta_i}{t_i}, which in turn depends on covariates via one or more smooth linear predictors with coefficients \eqn{\beta}{b}. The quadratic penalty terms penalize model complexity: \eqn{S_j}{S_j} is a known matrix and \eqn{\lambda_j}{v_j} an unknown smoothing parameter. Given smoothing parameters the penalized loss is readily minimized to estimate \eqn{\beta}{b}.

The smoothing parameters also have to be estimated. To this end, choose \eqn{k = 1,\ldots,m}{k=1...m} subsets \eqn{\alpha(k)\subset \{1,\ldots,n\}}{a(k) in \{1...n\}} and \eqn{\delta(k)\subset \{1,\ldots,n\}}{d(k) in \{1...n\}}. Usually \eqn{\delta(k)}{d(k)} is a subset of (or equal to) \eqn{\alpha(k)}{a(k)}. Let \eqn{\theta_i^{\alpha(k)}}{t_i^a(k)} denote the estimate of \eqn{\theta_i}{t_i} when the points indexed by \eqn{\alpha(k)}{a(k)} are omitted from fitting. Then the NCV criterion
\deqn{V = \sum_{k=1}^m \sum_{i \in \delta(k)} D(y_i,\theta_i^{\alpha(k)})}{V = sum_k sum_\{i in d(k)\}D(y_i,t_i^a(k)) }
is minimized w.r.t. the smoothing parameters, \eqn{\lambda_j}{v_j}. If \eqn{m=n}{m=n} and \eqn{\alpha(k)=\delta(k)=k}{a(k)=d(k)=k} then ordinary leave-one-out cross validation is recovered. This formulation covers many of the variants of cross validation reviewed in Arlot and Celisse (2010), for example.

Except for a quadratic loss, \eqn{V}{V} can not be computed exactly, but it can be computed to \eqn{O(n^{-2})}{O(n^-2)} accuracy (fixed basis size), by taking single Newton optimization steps from the full data \eqn{\beta}{b} estimates to the equivalent when each \eqn{\alpha(k)}{a(k)} is dropped. This is what \code{mgcv} does. The Newton steps require update of the full model Hessian to the equivalent when each datum is dropped. This can be achieved at \eqn{O(p^2)}{O(p^2)} cost, where \eqn{p}{p} is the dimension of \eqn{\beta}{b}. Hence, for example, the ordinary cross validation criterion is computable at the \eqn{O(np^2)}{O(np^2)} cost of estimating the model given smoothing parameters.

The NCV score computed in this way is optimized using a BFGS quasi-Newton method, adapted to the case in which smoothing parameters tending to infinity may cause indefiniteness.

For \code{bam(...,discrete=TRUE)} NCV can be used to estimate the smoothing parameters of the working penalized weighted linear model. This is typically very expensive relative to REML estimation, but the cost can be mitigated by basing the NCV criterion on a sample (usually random of the neighbourhoods). Note that in the case in which leave-out-neighbourhood CV is used to reduce the effects of autocorrelation, it is important to supply a neighbourhood for each observation, and to supply a \code{sample} element in \code{nei}. This ensures that the parameter covariance matrix is estimated correctly.     

}

\section{Spatial and temporal short range autocorrelation}{
A routine applied problem is that smoothing parameters tend to be underestimated in the presence of un-modelled short range autocorrelation, as the smooths try to fit the local excursions in the data caused by the local autocorrelation. Cross validation will tend to 'fit the noise' when there is autocorellation, since a model that fits the noise in the data correlated with an omitted datum, will also tend to closely fit the noise in the omitted datum, because of the correlation. That is autocorrelation works against the avoidance of overfit that cross validation seeks to achieve.

For short range autocorrelation the problems can be avoided, or at least mitigated, by predicting each datum when all the data in its `local' neighbourhood are omitted. The neighbourhoods being constructed in order that un-modelled correlation is minimized between the point of interest and points outside its neighbourhood. That is we set \eqn{m=n}{m=n}, \eqn{\delta(k)=k}{d(k)=k} and \eqn{\alpha(k) = {\tt nei}(k)}{a(k) = nei(k)}, where \code{nei(k)} are the indices of the neighbours of point \eqn{k}{k}. This approach has been known for a long time (e.g. Chu and Marron, 1991; Robert et al. 2017), but was previously rather too expensive for regular use for smoothing parameter estimation.}

\section{Specifying the neighbourhoods}{

The neighbourhood subsets \eqn{\alpha(k)}{a(k)} and \eqn{\delta(k)}{d(k)} have to be supplied to \code{\link{gam}}, and the \code{nei} argument does this. It is a list with the following arguments.
\itemize{
\item \code{a} is the vector of indices to be dropped for each neighbourhood.
\item \code{ma} gives the end of each neighbourhood. So \code{nei$a[(nei$ma[j-1]+1):nei$ma[j]]} gives the points dropped for the neighbourhood \code{j}: that is \eqn{\alpha(j)}{a(j)}.
\item \code{d} is the vector of indices of points to predict.
\item \code{md} gives the corresponding endpoints \code{md}. So \code{nei$d[(nei$md[j-1]+1):nei$md[j]]} indexes the points to predict for neighbourhood j: that is \eqn{\delta(j)}{d(j)}.
\item \code{sample} is an optional element used by \code{\link{bam}}. If it is a single number it gives the number of neighbourhoods to randomly sample to construct the NCV criterion. If it is a vector then it should contain the indices of the neighbourhoods to use.   
\item \code{jackknife} is an optional element used by \code{\link{gam}}. If supplied and \code{TRUE} then variance estimates are based on the raw Jackkife estimate, if \code{FALSE} then on the standard Bayesian results. If not supplied (usual) then an estimator accounting for the neighbourhood structure is used, that largely accounts for any correlation present within neighbourhoods. \code{jackknife} is ignored if NCV is being calculated for a model where another method is used for smoothing parameter selection.
}
If \code{nei==NULL} (or \code{a} or \code{ma} are missing) then leave-one-out cross validation is used. If \code{nei} is supplied but NCV is not selected as the smoothing parameter estimation method, then it is simply computed (but not optimized).
}

\section{Numerical issues}{
If a model is specified in which some coefficient values, \eqn{\beta}{b}, have non-finite likelihood then the NCV criterion computed with single Newton steps could also be non-finite. A simple fix replaces the NCV criterion with a quadratic approximation to the criterion around the full data fit. The quadratic approximation is always finite. This 'QNCV' is essential for some families, such as \code{\link{gevlss}}. QNCV is not needed for \code{\link{bam}} estimation.

Although the leading order cost of NCV is the same as REML or GCV, the actual cost is higher because the dominant operations costs are in matrix-vector, rather than matrix-matrix, operations, so BLAS speed ups are small. However multi-core computing is worthwhile for NCV. See the option \code{ncv.threads} in \code{\link{gam.control}}.
}




%- maybe also `usage' for other objects documented here.

\references{
Chu and Marron (1991) Comparison of two bandwidth selectors with dependent errors. The Annals of Statistics. 19, 1906-1918

Arlot, S. and A. Celisse (2010). A survey of cross-validation procedures for model selection. Statistics Surveys 4, 40-79

Roberts et al. (2017) Cross-validation strategies for data with temporal,
spatial, hierarchical, or phylogenetic structure. Ecography 40(8), 913-929.

Wood S.N. (2023) On Neighbourhood Cross Validation. \url{https://arxiv.org/abs/2404.16490}
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
require(mgcv)
nei.cor <- function(h,n) { ## construct nei structure
  nei <- list(md=1:n,d=1:n)
  nei$ma <- cumsum(c((h+1):(2*h+1),rep(2*h+1,n-2*h-2),(2*h+1):(h+1)))
  a0 <- rep(0,0); if (h>0) for (i in 1:h) a0 <- c(a0,1:(h+i))
  a1 <- n-a0[length(a0):1]+1
  nei$a <- c(a0,1:(2*h+1)+rep(0:(n-2*h-1),each=2*h+1),a1)
  nei
}
set.seed(1)
n <- 500;sig <- .6
x <- 0:(n-1)/(n-1)
f <- sin(4*pi*x)*exp(-x*2)*5/2
e <- rnorm(n,0,sig)
for (i in 2:n) e[i] <- 0.6*e[i-1] + e[i]
y <- f + e ## autocorrelated data
nei <- nei.cor(4,n) ## construct neighbourhoods to mitigate 
b0 <- gam(y~s(x,k=40)) ## GCV based fit
gc <- gam.control(ncv.threads=2)
b1 <- gam(y~s(x,k=40),method="NCV",nei=nei,control=gc)
## use "QNCV", which is identical here...
b2 <- gam(y~s(x,k=40),method="QNCV",nei=nei,control=gc)
## plot GCV and NCV based fits...
f <- f - mean(f)
par(mfrow=c(1,2))
plot(b0,rug=FALSE,scheme=1);lines(x,f,col=2)
plot(b1,rug=FALSE,scheme=1);lines(x,f,col=2)
}
\keyword{models} \keyword{regression}%-- one or more ..


