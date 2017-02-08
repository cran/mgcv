\name{gamlss.etamu}
\alias{gamlss.etamu}
\title{Transform derivatives wrt mu to derivatives wrt linear predictor}
\usage{
gamlss.etamu(l1, l2, l3 = NULL, l4 = NULL, ig1, g2, g3 = NULL,
  g4 = NULL, i2, i3 = NULL, i4 = NULL, deriv = 0)
}
\arguments{
\item{l1}{array of 1st order derivatives of log-likelihood wrt mu.}

\item{l2}{array of 2nd order derivatives of log-likelihood wrt mu.}

\item{l3}{array of 3rd order derivatives of log-likelihood wrt mu.}

\item{l4}{array of 4th order derivatives of log-likelihood wrt mu.}

\item{ig1}{reciprocal of the first derivative of the link function wrt the linear predictor.}

\item{g2}{array containing the 2nd order derivative of the link function wrt the linear predictor.}

\item{g3}{array containing the 3rd order derivative of the link function wrt the linear predictor.}

\item{g4}{array containing the 4th order derivative of the link function wrt the linear predictor.}

\item{i2}{two-dimensional index array, such that \code{l2[,i2[i,j]]} contains the partial w.r.t. params 
indexed by i,j with no restriction on the index values (except that they are in 1,...,ncol(l1)).}

\item{i3}{third-dimensional index array, such that \code{l3[,i3[i,j,k]]} contains the partial w.r.t. params 
indexed by i,j,k.}

\item{i4}{third-dimensional index array, such that \code{l4[,i4[i,j,k,l]]} contains the partial w.r.t. params 
indexed by i,j,k,l.}

\item{deriv}{if \code{deriv==0} only first and second order derivatives will be calculated. If
\code{deriv==1} the function goes up to 3rd order, and if \code{deriv==2} it provides
also 4th order derivatives.}
}
\value{
A list where the arrays \code{l1}, \code{l2}, \code{l3}, \code{l4} contain the derivatives (up
        to order four) of the log-likelihood wrt the linear predictor.
}
\description{Mainly intended for internal use in specifying location scale models.
Let \code{g(mu) = lp}, where \code{lp} is the linear predictor, and \code{g} is the link
         function. Assume that we have calculated the derivatives of the log-likelihood wrt \code{mu}.
         This function uses the chain rule to calculate the derivatives of the log-likelihood wrt
         \code{lp}. See \code{\link{trind.generator}} for array packing conventions. 
}
\seealso{\code{\link{trind.generator}}}

\author{
Simon N. Wood <simon.wood@r-project.org>.
}

