\name{gam.parser}
\alias{gam.parser}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Models fitting using penalized regression splines and 
GCV}
\description{ This is an internal function of package \code{mgcv}. It is a service routine for
 \code{gam} which splits a GAM formula into two parts: one representing the smooth part of the 
 model and the other representing the regular parametric part of the model (including any intercept term). 
 The processing of 
 the smooth terms in the model also returns the number of knots for each smooth term, and whether 
 or not the smooth term is to be treated as a pure regression spline (i.e. one without a penalty). 
 For further information on usuage see code for \code{gam}.
}
}



\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398


Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428


\url{http://www.ruwpa.st-and.ac.uk/simon.html}
}
\author{ Simon N. Wood \email{snw@st-and.ac.uk}}

\section{WARNING }{The code does not check for rank defficiency of the
model matrix that may result from including the same covariates in the smooth and parametric parts of the model -
      it will likely just fail instead!
 } 

\seealso{   \code{\link{gam}} \code{\link{mgcv}}}

\keyword{GAM, GCV, AIC, UBRE, penalized regression spline }%-- one or more ...

