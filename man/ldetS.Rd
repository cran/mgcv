\name{ldetS}
\alias{ldetS}
\title{Getting log generalized determinant of penalty matrices}
\usage{
ldetS(Sl, rho, fixed, np, root = FALSE,Stot=FALSE,repara = TRUE,
      nt = 1,deriv=2,sparse=FALSE)
}
\arguments{
\item{Sl}{the output of \code{Sl.setup}.}

\item{rho}{the log smoothing parameters.}

\item{fixed}{an array indicating whether the smoothing parameters are fixed (or free).}

\item{np}{number of coefficients.}

\item{root}{indicates whether or not to return the matrix square root, \code{E}, of the total penalty S_tot.}

\item{Stot}{indicates whether or not to return the total penalty S_tot.}

\item{repara}{if TRUE multi-term blocks will be re-parameterized using \code{gam.reparam}, and
a re-parameterization object supplied in the returned object.}

\item{nt}{number of parallel threads to use.}

\item{deriv}{order of derivative to use}
\item{sparse}{should \code{E} and/or \code{S} be sparse?}
}
\value{
A list containing: \itemize{
         \item{\code{ldetS}: the log-determinant of S. }
         \item{\code{ldetS1}: the gradient of the log-determinant of S. }
         \item{\code{ldetS2}: the Hessian of the log-determinant of S. }
         \item{\code{Sl}: with modified rS terms, if needed and rho added to each block }
         \item{\code{rp}:  a re-parameterization list. }
         \item{\code{E}:  a total penalty square root such that \code{t(E)\%*\%E = S_tot} (if \code{root==TRUE}). }
         \item{\code{S}:  the total penalty matrix (if \code{Stot==TRUE}).}

                         }
}
\description{
INTERNAL function calculating the log generalized determinant of penalty matrix S stored blockwise in an Sl list
             (which is the output of \code{Sl.setup}).
}
\author{
Simon N. Wood <simon.wood@r-project.org>.
}
