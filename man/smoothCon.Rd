\name{smoothCon}
\alias{smoothCon}
\alias{PredictMat}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction/Construction wrapper functions for GAM smooth terms}
\description{ Wrapper functions for construction of and prediction from smooth
terms in a GAM. The purpose of the wrappers is to allow user-transparant
re-parameterization of smooth terms, in order to allow identifiability
constraints to be absorbed into the parameterization of each term, if required.

}

\usage{
smoothCon(object,data,knots,absorb.cons=FALSE,scale.penalty=TRUE)
PredictMat(object,data)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{object}{ is a smooth specification object or a smooth object.}
\item{data}{A data frame containing the values of the (named) covariates at which the smooth term is to be 
evaluated.}
\item{knots}{An optional data frame supplying any knot locations to be
supplied for basis construction.}
\item{absorb.cons}{Set to \code{TRUE} in order to have identifiability
constraints absorbed into the basis.}
\item{scale.penalty}{should the penalty coefficient matrix be scaled to have
approximately the same `size' as the inner product of the terms model matrix
with itself? This can improve the performance of \code{\link{gamm}} fitting.}
}

\value{ From \code{smoothCon} a \code{smooth} object returned by the
appropriate \code{\link{smooth.construct}} method function. If constraints are
to be absorbed then the object will have an attributes \code{"qrc"} and
\code{"nCons"}, the qr decomposition of the constraint matrix (returned by
\code{\link{qr}}) and the number of constraints, respectively: these are  used in
the re-parameterization. 

For \code{predictMat} a matrix which will map the parameters associated with
the smooth to the vector of values of the smooth evaluated at the covariate
values given in \code{object}. 
 }


\details{ These wrapper functions exist to allow smooths specified using
\code{\link{smooth.construct}} and \code{\link{Predict.matrix}} method
functions to be re-parameterized so that identifiability constraints are no
longer required in fitting. This is done in a user transparent
manner, but is typically of no importance in use of GAMs. 

The parameterization used by \code{\link{gam}} can be controlled via
\code{\link{gam.control}}.
}


\references{ 

\url{http://www.maths.bath.ac.uk/~sw283/}
}

\author{Simon N. Wood \email{simon.wood@r-project.org}}

\seealso{ \code{\link{gam.control}},
\code{\link{smooth.construct}}, \code{\link{Predict.matrix}} }


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...

