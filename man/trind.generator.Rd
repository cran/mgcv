\name{trind.generator}
\alias{trind.generator}
\title{Generates index arrays for upper triangular storage}
\usage{
trind.generator(K = 2)
}
\arguments{
\item{K}{positive integer determining the size of the array.}
}
\value{
A list where the entries \code{i1} to \code{i4} are arrays in up to four dimensions, 
        containing K indexes along each dimension.
}
\description{
Generates index arrays for upper triangular storage up to order four. Useful when
             working with higher order derivatives, which generate symmetric arrays. 
             Mainly intended for internal use.
}
\details{
Suppose that \code{m=1} and you fill an array using code like
         \code{for(i in 1:K) for(j in i:K) for(k in j:K) for(l in k:K) 
          {a[,m] <- something; m <- m+1 }} and do this because actually the same 
          "something" would be stored for any permutation of the indices i,j,k,l.
          Clearly in storage we have the restriction l>=k>=j>=i, but for access we 
          want no restriction on the indices. \code{i4[i,j,k,l]} produces the 
          appropriate \code{m} for unrestricted indices. \code{i3} and {i2} do the same 
          for 3d and 2d arrays.
}
\examples{
library(mgcv)
A <- trind.generator(3)

# All permutations of c(1, 2, 3) point to the same index (5)
A$i3[1, 2, 3] 
A$i3[2, 1, 3]
A$i3[2, 3, 1]
A$i3[3, 1, 2]
A$i3[1, 3, 2]
}
\author{
Simon N. Wood <simon.wood@r-project.org>.
}

