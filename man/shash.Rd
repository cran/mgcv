\name{shash}
\alias{shash}
\title{Sinh-arcsinh location scale and shape model family}
\usage{
shash(link = list("identity", "logeb", "identity", "identity"), 
      b = 1e-2, phiPen = 1e-3)
}
\arguments{
\item{link}{vector of four characters indicating the link function for location, scale, skewness and kurtosis parameters.}
\item{b}{positive parameter of the logeb link function, see Details.}
\item{phiPen}{positive multiplier of a ridge penalty on kurtosis parameter. Do not touch it unless you know what you are doing, see Details.}
}
\value{
An object inheriting from class \code{general.family}.
}
\description{
The \code{shash} family implements the four-parameter sinh-arcsinh (shash) distribution of 
             Jones and Pewsey (2009). The location, scale, skewness and kurtosis of the density can depend 
             on additive smooth predictors. Useable only with gam, the linear predictors are specified 
             via a list of formulae. It is worth carefully considering whether the data are sufficient to support
	     estimation of such a flexible model before using it.
}
\details{
The density function of the shash family is 
         \deqn{p(y|\mu,\sigma,\epsilon,\delta)= C(z) \exp\{-S(z)^2/2\} \{2\pi(1+z^2)\}^{-1/2}/\sigma,}{p(y|\mu,\sigma,\epsilon,\delta)=C(z)exp{-S(z)^2/2}{2\pi(1+z^2)}^-0.5/\sigma}
         where \eqn{ C(z)=\{1+S(z)^2\}^{1/2} }{C(z)={1+S(z)^2}^0.5}, \eqn{ S(z)=\sinh\{\delta \sinh^{-1}(z)-\epsilon\} }{ S(z)=sinh{\delta asinh(z)-\epsilon} } and
         \eqn{z = (y - \mu)/(\sigma \delta)}. Here \eqn{\mu} and \eqn{\sigma > 0} control, respectively, location and 
         scale, \eqn{\epsilon} determines skewness, while \eqn{\delta > 0} controls tailweight.
         \code{shash} can model skewness to either side, depending on the sign of \eqn{\epsilon}. 
         Also, shash can have tails that are lighter (\eqn{\delta>1}) or heavier (\eqn{0<\delta<1}) that a normal.
         For fitting purposes, here we are using \eqn{\tau = \log(\sigma)}{\tau = log(\sigma)} and \eqn{\phi = \log(\delta)}{\phi = log(\delta)}.

The density is based on the expression given on the second line of section 4.1 and equation (2) of Jones and Pewsey (2009), and uses the simple reparameterization given in section 4.3. 
         
The link function used for \eqn{\tau} is logeb with is \eqn{\eta = \log \{\exp(\tau)-b\}}{\eta = log{exp(\tau)-b}} so that the inverse link is
\eqn{\tau = \log(\sigma) = \log\{\exp(\eta)+b\}}{\tau = log(\sigma) = log{exp(\eta)+b}}. The point is that we are don't allow \eqn{\sigma} to become smaller
than a small constant b. The likelihood includes a ridge penalty \eqn{- phiPen * \phi^2}, which shrinks \eqn{\phi} toward zero. When sufficient data is available the ridge penalty does not change the fit much, but it is useful to include it when fitting the model to small data sets, to avoid \eqn{\phi} diverging to +infinity (a problem already identified by Jones and Pewsey (2009)). 
}
\examples{

###############
# Shash dataset
###############
##  Simulate some data from shash
set.seed(847)
n <- 1000
x <- seq(-4, 4, length.out = n)

X <- cbind(1, x, x^2)
beta <- c(4, 1, 1)
mu <- X \%*\% beta 

sigma =  .5+0.4*(x+4)*.5            # Scale
eps = 2*sin(x)                      # Skewness
del = 1 + 0.2*cos(3*x)              # Kurtosis

dat <-  mu + (del*sigma)*sinh((1/del)*asinh(qnorm(runif(n))) + (eps/del))
dataf <- data.frame(cbind(dat, x))
names(dataf) <- c("y", "x")
plot(x, dat, xlab = "x", ylab = "y")

## Fit model
fit <- gam(list(y ~ s(x), # <- model for location 
                  ~ s(x),   # <- model for log-scale
                  ~ s(x),   # <- model for skewness
                  ~ s(x, k = 20)), # <- model for log-kurtosis
           data = dataf, 
           family = shash, # <- new family 
           optimizer = "efs")

## Plotting truth and estimates for each parameters of the density 
muE <- fit$fitted[ , 1]
sigE <- exp(fit$fitted[ , 2])
epsE <- fit$fitted[ , 3]
delE <- exp(fit$fitted[ , 4])

par(mfrow = c(2, 2))
plot(x, muE, type = 'l', ylab = expression(mu(x)), lwd = 2)
lines(x, mu, col = 2, lty = 2, lwd = 2)
legend("top", c("estimated", "truth"), col = 1:2, lty = 1:2, lwd = 2)

plot(x, sigE, type = 'l', ylab = expression(sigma(x)), lwd = 2)
lines(x, sigma, col = 2, lty = 2, lwd = 2)

plot(x, epsE, type = 'l', ylab = expression(epsilon(x)), lwd = 2)
lines(x, eps, col = 2, lty = 2, lwd = 2)

plot(x, delE, type = 'l', ylab = expression(delta(x)), lwd = 2)
lines(x, del, col = 2, lty = 2, lwd = 2)

## Plotting true and estimated conditional density
par(mfrow = c(1, 1))
plot(x, dat, pch = '.', col = "grey", ylab = "y", ylim = c(-35, 70))
for(qq in c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)){
  est <- fit$family$qf(p=qq, mu = fit$fitted)
  true <- mu + (del * sigma) * sinh((1/del) * asinh(qnorm(qq)) + (eps/del))
  lines(x, est, type = 'l', col = 1, lwd = 2)
  lines(x, true, type = 'l', col = 2, lwd = 2, lty = 2)
}
legend("topleft", c("estimated", "truth"), col = 1:2, lty = 1:2, lwd = 2)

#####################
## Motorcycle example
#####################

# Here shash is overkill, in fact the fit is not good, relative
# to what we would get with mgcv::gaulss
library(MASS)

b <- gam(list(accel~s(times, k=20, bs = "ad"), ~s(times, k = 10), ~1, ~1),
         data=mcycle, family=shash)

par(mfrow = c(1, 1))
xSeq <- data.frame(cbind("accel" = rep(0, 1e3),
                   "times" = seq(2, 58, length.out = 1e3)))
pred <- predict(b, newdata = xSeq)
plot(mcycle$times, mcycle$accel, ylim = c(-180, 100))
for(qq in c(0.1, 0.3, 0.5, 0.7, 0.9)){
  est <- b$family$qf(p=qq, mu = pred)
  lines(xSeq$times, est, type = 'l', col = 2)
}

plot(b, pages = 1, scale = FALSE)
}
\author{
Matteo Fasiolo <matteo.fasiolo@gmail.com> and Simon N. Wood.
}

\references{
Jones, M. and A. Pewsey (2009). Sinh-arcsinh distributions. Biometrika 96 (4), 761-780. \doi{10.1093/biomet/asp053}

Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\doi{10.1080/01621459.2016.1180986}
}

