\name{QT}
\alias{QT}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ QT factorisation of a matrix }
\description{ If \eqn{\bf A}{A} is an \eqn{r \times c}{r by c} matrix, 
where \eqn{r<c}{r<c}, then it can be factorised: 
\eqn{ {\bf AQ}=[{\bf 0},{\bf T}] }{ AQ=[0,T] }, 
where \eqn{\bf Q}{Q} is an orthogonal matrix and \eqn{\bf T}{T}
is reverse lower triangular. The first \eqn{c-r}{c-r} columns of \eqn{\bf Q}{Q} 
from an orthogonal basis for the null space of \eqn{\bf A}{A}. \eqn{\bf Q}{Q} is 
returned as a sequence of \eqn{r}{r} Householder rotations (from the right), as 
described below.
}
\usage{
QT(A)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{A}{ is an \R matrix having more columns than rows }
}
\details{ This function is primarily useful for providing the null space of the 
linear constraint matrix \eqn{\bf C}{C}, from the linear constraints 
\eqn{ {\bf Cp} = {\bf 0} }{ Cp=0 } applying to a genralized ridge regression problem to be
solved by 
\code{mgcv()}. The format of the null space as provided by the returned matrix is 
as expected by \code{mgcv()}. 

   
}
\value{

The returned matrix \eqn{\bf R}{R}, say, is of the same dimension as the input 
matrix \eqn{\bf A}{A}. Each of its rows contains a vector, \eqn{ {\bf u}_i }{ u_i }, 
defining one Householder rotation, \eqn{ {\bf H}_i =({\bf I} - {\bf u}_i {\bf 
u}_i^\prime) }{ H_i = (I - u_i u_i') }. The orthogonal matrix \eqn{\bf Q}{Q} is 
defined by: \eqn{ {\bf Q}={\bf H}_1 {\bf H}_2 \ldots {\bf H}_r }{ Q=H_1 H_2 ... 
H_r }.
 
}
\references{Gill, Murray and Wright (1981) Practical Optimization, Academic %%@
Press}

\author{ Simon N. Wood \email{snw@st-and.ac.uk}}


\examples{
y<-rnorm(12,1,1)
A<-matrix(y,2,6)
B<-QT(A)
}
\keyword{ QT orthogonal factorisation }%-- one or more ...
