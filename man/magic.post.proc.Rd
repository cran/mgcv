\name{magic.post.proc}
\alias{magic.post.proc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Auxilliary information from magic fit}
\description{Obtains parameter covariance matrix, estimated degrees of 
freedom for each parameter and leading diagonal of influence/hat matrix, 
for a penalized regression estimated by \code{magic}.
}
\usage{
magic.post.proc(X,object,w)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{X}{ is the model matrix.}
\item{object}{is the list returned by \code{magic} after fitting the 
model with model matrix \code{X}.}
\item{w}{is the weight vector used in fitting, or the weight matrix used 
in fitting (i.e. supplied to 
\code{magic}, if one was.)}
}
\details{ \code{object} contains \code{rV} (\eqn{ {\bf V}}{V}, say), and 
\code{scale} (\eqn{ \phi}{s}, say) which can be 
used to obtain the require quantities as follows. The covariance matrix of 
the parameters is \eqn{ {\bf VV}^\prime \phi}{VV's}. The vector of 
estimated degrees of freedom for each parameter is the leading diagonal of 
\eqn{ {\bf VV}^\prime {\bf X}^\prime {\bf W}^\prime {\bf W}{\bf X}}{ VV'X'W'WX} 
where \eqn{\bf{W}}{W} is either the 
weight matrix \code{w} or the matrix \code{diag(w)}. The 
hat/influence  matrix is given by 
\eqn{ {\bf WX}{\bf VV}^\prime {\bf X}^\prime {\bf W}^\prime }{ WXVV'X'W'} 
.
}
\value{ A list with three items:
\item{Vb}{the covariance matrix of the model parameters.}
\item{hat}{the leading diagonal of the hat (influence) matrix.}
\item{edf}{the array giving the estimated degrees of freedom associated 
with each parameter.}
}
 
\seealso{\code{\link{magic}}}

\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









