\name{GAMsetup}
\alias{GAMsetup}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Set up GAM using penalized cubic regression splines} 
\description{ Sets up design matrix \eqn{\bf X}{X}, penalty matrices \eqn{ {\bf
S}_i }{S_i} and constraint matrix null space \eqn{\bf Z}{Z} for a GAM defined in terms of 
penalized regression splines, as well as returning the locations of the knots of these 
regression splines \code{xp[][]}. The output is such that the model can be fitted and 
smoothing parameters estimated by the method of Wood (2000) as implemented in routine
\code{mgcv()}.

}
\usage{
GAMsetup(G)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  The function takes a single argument \code{G}, but this is a list containing several 
elements:
\item{m}{the number of smooth terms in the model}
\item{n}{the number of data to be modelled}
\item{nsdf}{the number of user supplied columns of the design matrix for any parametric 
model parts}
\item{df}{an array of \code{G$m} integers specifying the maximum d.f. for each spline 
term.}
\item{x}{an array of \code{G$n} element arrays of data and (optionally) design matrix 
columns. The first \code{G$nsdf} elements of \code{G$x} should contain the elements of 
the columns of the design matrix corresponding to the parametric part of the model. The 
remaining \code{G$m} elements of \code{G$x} are the values of the covariates that are 
arguments of the spline terms. Note that the user should not supply a column of 1's in 
the design matrix, corresponding to the mean response as this will be added 
automatically.}
}
\value{
  A list \code{H}, containing the elements of \code{G} (the input list) plus the 
following:  
\item{X}{the full design matrix.}
\item{S}{an array of matrices containing the coefficients of the penalties. These are 
stored in a compact form, so that \code{H$S[i]} is the smallest square matrix containing 
submatrix  containing all the non-zero elements of \eqn{ {\bf S}_i}{S_i}
the ith penalty 
matrix. Element 0,0 of \code{H$S[i]} is element \code{off[i],off[i]} of \eqn{ {\bf 
S}_i }{S_i}, element 0,1 of \code{H$S[i]} is element \code{off[i],off[i]+1} of \eqn{ {\bf 
S}_i }{S_i}, and so on.}
\item{off}{is an array of offsets, used to facilitate efficient storage of the penalty 
matrices and to indicate where in the overall parameter vector the parameters of the ith 
spline reside (e.g. first parameter of ith spline is at \code{p[off[i]]}).}
\item{Z}{null space of constraints stored as series of Householder rotations.}  
\item{xp}{matrix whose rows contain the covariate values corresponding to the 
parameters  of each spline - the splines are parameterized using their \eqn{y}{y}- 
values at a series of \eqn{x}{x} values - these vectors contain those \eqn{x}{x} 
values!}
}
\references{ Wood, S.N. (2000) Modelling and smoothing parameter estimation with multiple
quadratic penalties" JRSSB 62(2) }
\author{Simon N. Wood \email{snw@st-and.ac.uk}}

\seealso{ \code{\link{mgcv}} }

\examples{
    # This example modified from routine SANtest()

    n<-100 # number of observations to simulate
    x <- runif(4 * n, 0, 1) # simulate covariates
    x <- array(x, dim = c(4, n)) # put into array for passing to GAMsetup
    pi <- asin(1) * 2  # begin simulating some data
    y <- 2 * sin(pi * x[1, ])
    y <- y + exp(2 * x[2, ]) - 3.75887
    y <- y + 0.2 * x[3, ]^11 * (10 * (1 - x[3, ]))^6 + 10 * (10 * 
        x[3, ])^3 * (1 - x[3, ])^10 - 1.396
    sig2<- -1    # set magnitude of variance 
    e <- rnorm(n, 0, sqrt(abs(sig2)))
    y <- y + e          # simulated data
    w <- matrix(1, n, 1) # weight matrix
    par(mfrow = c(2, 2)) # scatter plots of simulated data
    plot(x[1, ], y)
    plot(x[2, ], y)
    plot(x[3, ], y)
    plot(x[4, ], y)
    G <- list(m = 4, n = n, nsdf = 0, df = c(15, 15, 15, 15), 
        x = x) # creat list for passing to GAMsetup
    H <- GAMsetup(G)
    H$y <- y    # add data to H
    H$sig2 <- sig2  # add variance (signalling GCV use in this case) to H
    H$w <- w       # add weights to H
    H <- mgcv(H)  # select smoothing parameters and fit model
    
}
\keyword{GAM }%-- one or more ...
