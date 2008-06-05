\name{smoothCon}
\alias{smoothCon}
\alias{PredictMat}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction/Construction wrapper functions for GAM smooth terms}
\description{ Wrapper functions for construction of and prediction from smooth
terms in a GAM. The purpose of the wrappers is to allow user-transparant
re-parameterization of smooth terms, in order to allow identifiability
constraints to be absorbed into the parameterization of each term, if required.
The routine also handles `by' variables and construction of identifiability constraints automatically, 
although this behaviour can be over-ridden.

}

\usage{
smoothCon(object,data,knots,absorb.cons=FALSE,scale.penalty=TRUE,n=nrow(data),dataX=NULL)
PredictMat(object,data,n=nrow(data))
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{object}{ is a smooth specification object or a smooth object.}
\item{data}{A data frame, model frame or list containing the values of the 
(named) covariates at which the smooth term is to be 
evaluated. If it's a list then \code{n} must be supplied.}
\item{knots}{An optional data frame supplying any knot locations to be
supplied for basis construction.}
\item{absorb.cons}{Set to \code{TRUE} in order to have identifiability
constraints absorbed into the basis.}
\item{scale.penalty}{should the penalty coefficient matrix be scaled to have
approximately the same `size' as the inner product of the terms model matrix
with itself? This can improve the performance of \code{\link{gamm}} fitting.}
\item{n}{number of values for each covariate, or if a covariate is a matrix, 
the number of rows in that matrix: must be supplied explicitly if \code{data} is a list. }
\item{dataX}{Sometimes the basis should be set up using data in \code{data}, but the model matrix
should be constructed with another set of data provided in \code{dataX} --- \code{n} is assumed to 
be the same for both. Facilitates smooth id's.}
}

\value{ From \code{smoothCon} a list of \code{smooth} objects returned by the
appropriate \code{\link{smooth.construct}} method function. If constraints are
to be absorbed then the objects will have  attributes \code{"qrc"} and
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
manner, but is typically of no importance in use of GAMs. The routine's 
also handle \code{by} variables and will create default identifiability 
constraints.

If a user defined smooth constructor handles \code{by} variables itself, then its 
returned smooth object should contain an object \code{by.done}. If this does not exist 
then \code{smoothCon} will use the default code. Similarly if a user defined \code{Predict.matrix} 
method handles \code{by} variables internally then the returned matrix should have a 
\code{"by.done"} attribute.

Default centering constraints, that terms should sum to zero over the covariates, are produced unless 
the smooth constructor includes a matrix \code{C} of constraints. To have no constraints (in which case 
you had better have a full rank penalty!) the matrix \code{C} should have no rows.

\code{smoothCon} returns a list of smooths because factor \code{by} variables result in multiple copies 
of a smooth, each multiplied by the dummy variable associated with one factor level. \code{smoothCon} modifies 
the smooth object labels in the presence of \code{by} variables, to ensure that they are unique, it also stores 
the level of a by variable factor associated with a smooth, for later use by \code{PredMat}.

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

