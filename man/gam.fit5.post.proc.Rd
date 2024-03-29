
\name{gam.fit5.post.proc}
\alias{gam.fit5.post.proc}
\title{Post-processing output of gam.fit5}
\usage{
gam.fit5.post.proc(object, Sl, L, lsp0, S, off, gamma)
}
\arguments{
\item{object}{output of \code{gam.fit5}.}

\item{Sl}{penalty object, output of \code{Sl.setup}.}

\item{L}{matrix mapping the working smoothing parameters.}

\item{lsp0}{log smoothing parameters.}

\item{S}{penalty matrix.}

\item{off}{vector of offsets.}

\item{gamma}{parameter for increasing model smoothness in fitting.}
}
\value{
A list containing: \itemize{
         \item{\code{R}: unpivoted Choleski of estimated expected hessian of log-likelihood. }
         \item{\code{Vb}: the Bayesian covariance matrix of the model parameters. }
         \item{\code{Ve}: "frequentist" alternative to \code{Vb}. }
         \item{\code{Vc}: corrected covariance matrix. }
         \item{\code{F}: matrix of effective degrees of freedom (EDF). }
         \item{\code{edf}:  \code{diag(F)}. }
         \item{\code{edf2}:  \code{diag(2F-FF)}. } 
                         }
}
\description{
INTERNAL function for post-processing the output of \code{gam.fit5}.
}
\author{
Simon N. Wood <simon.wood@r-project.org>.
}
