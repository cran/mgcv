\name{blas.thread.test}
\alias{blas.thread.test}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{BLAS thread safety}
\description{Most BLAS implementations are thread safe, but some versions of OpenBLAS, for example, are not. This routine is a diagnostic helper function, which you will never need if you don't set \code{nthreads>1}, and even then are unlikely to need.   
}
\usage{
blas.thread.test(n=1000,nt=4)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{n}{Number of iterations to run of parallel BLAS calling code.}
 \item{nt}{Number of parallel threads to use}
}

\details{While single threaded OpenBLAS 0.2.20 was thread safe, versions 0.3.0-0.3.6 are not, and from version 0.3.7 thread safety of the single threaded OpenBLAS requires making it with the option \code{USE_LOCKING=1}. The reference BLAS is thread safe, as are MKL and ATLAS. This routine repeatedly calls the BLAS from multi-threaded code and is sufficient to detect the problem in single threaded OpenBLAS 0.3.x.

A multi-threaded BLAS is often no faster than a single-threaded BLAS, while judicious use of threading in the code calling the BLAS can still deliver a modest speed improvement. For this reason it is often better to use a single threaded BLAS and the code{nthreads} options to \code{\link{bam}} or \code{\link{gam}}. For \code{bam(...,discrete=TRUE)} using several threads can be a substantial benefit, especially with the reference BLAS.

The MKL BLAS is mutlithreaded by default. Under linux setting environment variable \code{MKL_NUM_THREADS=1} before starting R gives single threaded operation. 

}
 



\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









