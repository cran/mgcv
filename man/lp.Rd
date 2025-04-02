\name{lp}
\alias{lp}
\alias{feasible}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Basic linear programming}

\description{Solves, for \eqn{{\bf x}}{x}, linear programming problems in the standard form
\deqn{\min {\bf c}^T {\bf x} \text{ s.t. } {\bf Ax}={\bf b},~~{\bf b}\ge {\bf 0}}{min c'x s.t. Ax=b, x>=0}
or in the form
\deqn{\min {\bf c}^T {\bf x} \text{ s.t. } {\bf Ax}\ge {\bf b},~~{\bf Cx}= {\bf d}}{min c'x s.t. Ax>=b, Cx=d}
A lightweight implementation of the simplex method, designed for problems of modest size, not large scale problems requiring sparse methods.  

\code{feasible} finds starting values meeting the constraints, which is also useful for finding feasible initial values for quadratic programming.  
}
\usage{
lp(c,A,b,C=NULL,d=NULL,Bi=NULL,maxit=max(1000, nrow(A) * 10), phase1 = FALSE)
feasible(A,b,C=NULL,d=NULL,maxit = max(1000, nrow(A) * 10))
}
%- maybe also `usage' for other objects documented here.

\arguments{ 
\item{c}{The vector defining the linear program objective function \eqn{{\bf c}^T{\bf x}}{c'x}.} 

\item{A}{The constraint matrix.} 

\item{b}{The vector defining the r.h.s. of the constraint involving \code{A}.} 

\item{C}{\code{NULL} to signal a problem in standard form. Otherwise the matrix defining the equality constraints.}

\item{d}{\code{NULL} to signal a problem in standard form. Otherwise the vector defining the r.h.s. of the equality constraint.}

\item{Bi}{For a problem in standard form, index of the elements of \eqn{\bf x}{x} initially non-zero. \eqn{{\bf A}[,Bi]{\bf x}^* = {\bf b}}{A[,Bi]x*=b} must have a unique solution with \eqn{{\bf x}^*\ge {\bf 0}}{x*>=0}. \code{NULL} to have this set found automatically by solution of a phase 1 linear program.}

\item{maxit}{Maximum number of iterations of the simplex method to allow.}

\item{phase1}{signals that function is being called to solve a phase 1 problem.}

}

\value{A vector. Either the solution of the problem (\code{lp}), or a feaible initial vector (\code{feasible}). Produces an error if there is no solution or \code{maxit} is exceeded.
}


\details{This code was written to provide a lightweight method for finding feasible initial coefficients for shape constrained splines, but is suitable for solving general linear programming problems of modest size. Given that it uses dense matrix computations, it is not suitable for large scale problems where it is important to exploit sparcity. Neither is it suitable for the sort of linear programming problems arising from integer programs, where degeneracy may be a substantial problem.   

The function uses the simplex method (see e.g. Chapter 13 of Nocedal and Wright, 2006). Computational efficiency is ensured by QR decomposing \eqn{{\bf A}[,Bi]}{A[,Bi]} and then efficiently updating the decomposition, every time the active set \code{Bi} changes, using Givens based updating schemes broadly similar to those given in Golub and van Loan (2013) 5.1.3 p240 and 6.5.3 p337. Given the QR factorization solution of systems involving \eqn{{\bf A}[,Bi]}{A[,Bi]} is efficient. }


\references{
Golub G.H. and C.F. van Loan (2013) Matrix Computations (4th edition) Johns Hopkins

Nocedal J. and S. Wright (2006) Numerical Optimization (2nd edition) Springer 
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\section{WARNINGS }{
Not designed for large scale problems requiring sparse methods, nor for problems where significant degeneracy is expected. 
} 

\seealso{\code{\link{pcls}}
 }

\examples{
library(mgcv)
## very simple linear program...
c <- c(-4,-2,0,0)
A <- matrix(c(1,2,1,.5,1,0,0,1),2,4)
b <- c(5,8)
x <- lp(c,A,b,c(3,4));sum(c*x);x ## Bi given

x <- lp(c,A,b);sum(c*x);x ## Bi found automatically

## equivalent formulation Ax>b
A <- -matrix(c(1,2,1,.5),2,2)
b <- -c(5,8)
C <- matrix(0,0,2); d <- numeric(0)
c <- c(-4,-2)
x <- lp(c,A,b,C=C,d=d);sum(c*x);x

## equivalent formulation Ax>b, Cx=d
c <- c(-4,-2,0)
A <- matrix(c(0,-2,0,-.5,1,0),2,3)
C <- matrix(c(1,1,1),1,3); d <- 5
b <- c(0,-8)
x <- lp(c,A,b,C=C,d=d);sum(c*x);x
}

\keyword{linear program} \keyword{shape constraint}






