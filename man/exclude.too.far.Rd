\name{exclude.too.far}
\alias{exclude.too.far}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Exclude prediction grid points too far from data}
\description{ Takes two arrays defining the nodes of a grid over a 2D covariate space and two arrays 
defining the location of data in that space, and returns a logical vector with elements \code{TRUE} if 
corresponding node is too far from data and \code{FALSE} otherwise. Basically a service routine for 
\code{vis.gam}.
}
\usage{
exclude.too.far(g1,g2,d1,d2,dist)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{g1}{co-ordinates of grid relative to first axis.}
\item{g2}{co-ordinates of grid relative to second axis.}
\item{d1}{co-ordinates of data relative to first axis.}
\item{d2}{co-ordinates of data relative to second axis.}
\item{dist}{how far away counts as two far. Grid and data are first scaled so that the grid lies exactly 
in the unit square, and \code{dist} is a distance within this unit square.} 
}
\details{ Linear scalings of the axes are first determined so that the grid defined by the nodes in 
\code{g1} and \code{g2} lies exactly in the unit square (i.e. on [0,1] by [0,1]). These scalings are 
applied to \code{g1}, \code{g2}, \code{d1} and \code{d2}. The minimum Euclidean 
distance from each node to a datum is then determined and if it is greater than \code{dist} the 
corresponding entry in the returned array is set to \code{TRUE} (otherwise to \code{FALSE}). The 
calculations are performed by looping through the nodes (although they are vectorized with respect to the 
data): unless there are very few data this is not overly ineffeicient, and avoids using excessive amounts of 
memory. 
}

\value{An logical array with \code{TRUE} indicating a node in the grid defined by \code{g1}, \code{g2} that 
is too far from any datum to be included. 
}

\references{
\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk}}

\section{WARNING }{ Routine is vectorized w.r.t. the data and not the grid: it will be very inefficient if 
\code{d1} and \code{d2} are very short while \code{g1} and \code{g2} are very long.
} 

\seealso{ \code{\link{vis.gam}} }

\examples{
library(mgcv)
x<-rnorm(100);y<-rnorm(100) # some "data"
n<-20 # generate a grid....
mx<-seq(min(x),max(x),length=n)
my<-seq(min(y),max(y),length=n)
gx<-rep(mx,n);gy<-rep(my,rep(n,n))
tf<-exclude.too.far(gx,gy,x,y,0.1)
plot(gx[!tf],gy[!tf],pch=".");points(x,y,col=2)
}
\keyword{hplot}%-- one or more ...






