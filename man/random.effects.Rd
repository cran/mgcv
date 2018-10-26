\name{random.effects}
\alias{random.effects}
\title{Random effects in GAMs}
\description{ The smooth components of GAMs can be viewed as random effects for estimation purposes. This means that more conventional 
random effects terms can be incorporated into GAMs in two ways. The first method converts all the smooths into fixed and random components 
suitable for estimation by standard mixed modelling software. Once the GAM is in this form then conventional random effects are easily added,
and the whole model is estimated as a general mixed model. \code{\link{gamm}} and \code{gamm4} from the \code{gamm4} package operate in this way.

The second method represents the conventional random effects in a GAM in the same way that the smooths are represented --- as penalized 
regression terms. This method can be used with \code{\link{gam}} by making use of \code{s(...,bs="re")} terms in a model: see 
\code{\link{smooth.construct.re.smooth.spec}}, for full details. The basic idea is that, e.g., \code{s(x,z,g,bs="re")} generates an i.i.d. Gaussian 
random effect with model matrix given by \code{model.matrix(~x:z:g-1)} --- in principle such terms can take any number of arguments. This simple 
approach is sufficient for implementing a wide range of commonly used random effect structures. For example if \code{g} is a factor then 
\code{s(g,bs="re")} produces a random coefficient for each level of \code{g}, with the random coefficients all modelled as i.i.d. normal. If \code{g} is a factor and \code{x} is numeric, then \code{s(x,g,bs="re")} produces an i.i.d. normal random slope relating the response to \code{x} for each level of \code{g}. If \code{h} is another factor then \code{s(h,g,bs="re")} produces the usual i.i.d. normal \code{g} - \code{h} interaction. Note that a rather useful approximate test for zero random effect is also implemented for such terms based on Wood (2013). If the precision
matrix is known to within a multiplicative constant, then this can be supplied via the \code{xt} argument of \code{s}. See \link{smooth.construct.re.smooth.spec} for details and example. 

Alternatively, but less straightforwardly, the \code{paraPen} argument to \code{\link{gam}} can be used: 
see \code{\link{gam.models}}. If smoothing parameter estimation is by ML or REML (e.g. \code{gam(...,method="REML")}) then this approach is 
a completely conventional likelihood based treatment of random effects.

\code{gam} can be slow for fitting models with large numbers of random effects, because it does not exploit the sparsity that is often a feature
of parametric random effects. It can not be used for models with more coefficients than data. However \code{gam} is often faster and more reliable 
than \code{gamm} or \code{gamm4}, when the number of random effects is modest.

To facilitate the use of random effects with \code{gam}, \code{\link{gam.vcomp}} is a utility routine for converting 
smoothing parameters to variance components. It also provides confidence intervals, if smoothness estimation is by ML or REML. 

Note that treating random effects as smooths does not remove the usual problems associated with testing variance components for equality to zero: see \code{\link{summary.gam}} and \code{\link{anova.gam}}.   

}

\seealso{\code{\link{gam.vcomp}}, \code{\link{gam.models}}, \code{\link{smooth.terms}}, 
\code{\link{smooth.construct.re.smooth.spec}},
\code{\link{gamm}}}

\author{
Simon Wood <simon.wood@r-project.org>
}
\references{
Wood, S.N. (2013) A simple test for random effects in regression models. Biometrika 100:1005-1010

Wood, S.N. (2011) Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

Wood, S.N. (2008) Fast stable direct fitting and smoothness
 selection for generalized additive models. Journal of the Royal
 Statistical Society (B) 70(3):495-518

Wood, S.N. (2006) Low rank scale invariant tensor product smooths for
generalized additive mixed models. Biometrics 62(4):1025-1036

}

\examples{
## see also examples for gam.models, gam.vcomp, gamm
## and smooth.construct.re.smooth.spec

## simple comparison of lme and gam
require(mgcv)
require(nlme)
b0 <- lme(travel~1,data=Rail,~1|Rail,method="REML") 

b <- gam(travel~s(Rail,bs="re"),data=Rail,method="REML")

intervals(b0)
gam.vcomp(b)
anova(b)
plot(b)

## simulate example...
dat <- gamSim(1,n=400,scale=2) ## simulate 4 term additive truth

fac <- sample(1:20,400,replace=TRUE)
b <- rnorm(20)*.5
dat$y <- dat$y + b[fac]
dat$fac <- as.factor(fac)

rm1 <- gam(y ~ s(fac,bs="re")+s(x0)+s(x1)+s(x2)+s(x3),data=dat,method="ML")
gam.vcomp(rm1)

fv0 <- predict(rm1,exclude="s(fac)") ## predictions setting r.e. to 0
fv1 <- predict(rm1) ## predictions setting r.e. to predicted values

}
\keyword{regression}