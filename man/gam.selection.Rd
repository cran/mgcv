\name{gam.selection}
\alias{gam.selection}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized Additive Model Selection}
\description{ This page is intended to provide some more information on
  how to select GAMs. Given a model structure specified by a gam model formula,
\code{gam()} attempts to find the appropriate smoothness for each applicable model 
term using Generalized Cross Validation (GCV) or an Un-Biased Risk Estimator (UBRE),
the latter being used in cases in which the scale parameter is assumed known. GCV and 
UBRE are covered in Craven and Wahba (1979) and Wahba (1990). The method used alternates 
grid searches for the correct overall level of smoothness for the whole model, given the 
relative smoothness of terms, with Newton/Steepest descent updates of the relative smoothness 
of terms, given the overall amount of smoothness. 

Automatic smoothness selection is unlikely to be successful with few data, particularly with
multiple terms to be selected. The methods can also fail to find the real minimum of the GCV/UBRE
score if the model contains many smooth terms that should really be completely smooth, or close to it (e.g. a straight line 
for a default 1-d smooth). The problem is that in this circumstance the optimal overall smoothness given the relative
smoothness of terms may make all terms completely smooth - but this will tend to move the smoothing parameters 
to a location where the GCV/UBRE score is nearly completely flat with respect to the smoothing parameters so that Newton and steepest 
descent are both ineffective. These problems can usually be overcome by replacing some completely smooth terms with purely
parametric model terms. 

A good example of where smoothing parameter selection can ``fail'', but in an unimportant manner is provided by the 
\code{rock.gam} example in Venables and Ripley. In this case 3 smoothing parameters are to estimated from 48 data, which is 
probably over-ambitious. \code{gam} will estimate either 1.4 or 1 degrees of freedom for the smooth of \code{shape}, depending on 
the exact details of model specification (e.g. k value for each \code{s()} term). The lower GCV score is really at 1.4 (and if the other
2 terms are replaced by straight lines this estimate is always returned), but the \code{shape} term is in no way significant and the 
lowest GCV score is obtained by removing it altogether. The problem here is that the GCV score contains very little information on the optimal 
degrees of freedom to associate with a term that GCV would suggest should really be dropped. 


In general the most logically consistent method to use for deciding which terms to include in the model is to compare GCV/UBRE scores
for models with and without the term. More generally the score for the model with a smooth term can be compared to the score for the model with 
the smooth term replaced by appropriate parametric terms. Candidates for removal can be identified by reference to the approximate p-values provided by \code{summary.gam}. Candidates for replacement by parametric terms are smooth terms with estimated degrees of freedom close to their 
minimum possible.

}

\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\references{
Craven and Wahba (1979) Smoothing Noisy Data with Spline Functions. Numer. Math. 31:377-403

Venables and Ripley (1999) Modern Applied Statistics with S-PLUS

Wahba (1990) Spline Models of Observational Data. SIAM.

Wood (2000) Modelling and Smoothing Parameter Estimation 
with Multiple 
   Quadratic Penalties. JRSSB 62(2):413-428

}

\keyword{models} \keyword{regression}%-- one or more ..


