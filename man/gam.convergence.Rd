\name{gam.convergence}
\alias{gam.convergence}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM convergence issues.}
\description{ Sometimes the GAM estimation process fails to converge. There
are two possible reasons for this (assuming that the obvious possibility of
insufficient iterations being allowed by \code{gam.control} has been
excluded). 
Firstly, a GAM is just a (penalized) GLM
and the IRLS scheme used to estimate GLMs is not guaranteed to
converge. Hence non convergence of a GAM may relate to a lack of stability in
the basic IRLS scheme. Secondly,  nesting GAM smoothing parameter estimation 
within the IRLS iteration loop (the mgcv default) can also make convergence
more difficult to achieve. Of course, these two scenarios are not really
completely distinct, but thinking about convergence in these terms does lead
to a useful practical approach for achieving convergence. 

Basically, the first step is to establish whether the IRLS iterations
are capable of converging. To do this fit the problematic GAM with all smooth
terms specified with \code{fx=TRUE} so that the smoothing parameters are all
fixed at zero. If this `largest' model can converge then it is the smoothing
parameter estimation that is causing the problem and convergence can usaully
be achieved by setting the \code{gam.control} option
\code{spIterType="outer"}. This iterates the IRLS scheme to convergence for
each trial set of smoothing parameters, and is usually much slower than the
default \code{spIterType="perf"}, but it avoids the smoothing parameter
estimation de-stabilizing the IRLS iterations. Note that this approach does
not in anyway change the model being fitted.

If, on the other hand, the model with all smoothing parameters set to zero does
not converge, then the problem lies with the IRLS process itself. Having tried
increasing \code{maxit} in \code{gam.control}, there are several other
possibilities for stabilizing the iteration. All are based on the fact that
too much model flexibility is almost always a major factor in
non-convergence. It is possible to try (i) setting lower bounds on the
smoothing parameters using the \code{min.sp} argument of \code{gam}: this may
or may not change the model being fitted; (ii)
reducing the flexibility of the model by reducing the basis dimensions
\code{k} in the specification of \code{s} and \code{te} model terms: this
obviously changes the model being fitted somewhat; (iii)
introduce a small regularization term into the fitting via the \code{irls.reg}
argument of \code{gam.control}: this option obviously changes the nature of
the fit somewhat, since parameter estimates are pulled towards zero by doing
this. 
Once the basic IRLS process has been stabilized, you can try calling
\code{gam} using the stabilized model and its
options. If the default fitting fails again, then retry using the
\code{spIterType="outer"} option discussed previously. 

In most cases model fitting can be made to converge in this way, and the
stabilization often has little or no real affect on the estimates. However
care and judgement are needed, since some of these techniques have the
potential to have a large impact on estimates. 

Often, of course, a major contributer to fitting difficulties is that the
model is a very poor description of the data.  
 

}

%- maybe also `usage' for other objects documented here.


\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
