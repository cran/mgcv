\name{trind.generator}
\alias{trind.generator}
\title{Generates index arrays for upper triangular storage}
\usage{
trind.generator(K = 2, ifunc=FALSE, reverse= !ifunc)
}
\arguments{
\item{K}{positive integer determining the size of the array.}
\item{ifunc}{if \code{TRUE} index functions are returned in place of index arrays.}
\item{reverse}{should the reverse indices be computed? Probably not if \code{ifunc==TRUE}.}
}
\value{
A list where the entries \code{i1} to \code{i4} are arrays in up to four dimensions, 
        containing K indexes along each dimension. If \code{ifunc==TRUE} index functions
	are returned in place of index arrays. If \code{reverse==TRUE} reverse indices
	\code{i1r} to \code{i4r} are returned (always as arrays).
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
          appropriate \code{m} for unrestricted indices. \code{i3} and \code{i2} do the same 
          for 3d and 2d arrays. If \code{ifunc==TRUE} then \code{i2}, \code{i3} and \code{i4}
	  are functions, so \code{i4(i,j,k,l)} returns appropriate \code{m}. For high \code{K}
	  the function versions save storage, but are slower.

If computed, the reverse indices pick out the unique elements of a symmetric array stored redundantly.
The indices refer to the location of the elements when the redundant array is accessed as its underlying
vector. For example the reverse indices for a 3 by 3 symmetric matrix are 1,2,3,5,6,9.
}
\examples{
library(mgcv)
A <- trind.generator(3,reverse=TRUE)

# All permutations of c(1, 2, 3) point to the same index (5)
A$i3[1, 2, 3] 
A$i3[2, 1, 3]
A$i3[2, 3, 1]
A$i3[3, 1, 2]
A$i3[1, 3, 2]

## use reverse indices to pick out unique elements
## just for illustration...
A$i2;A$i2[A$i2r]
A$i3[A$i3r]


## same again using function indices...
A <- trind.generator(3,ifunc=TRUE)
A$i3(1, 2, 3) 
A$i3(2, 1, 3)
A$i3(2, 3, 1)
A$i3(3, 1, 2)
A$i3(1, 3, 2)
}
\author{
Simon N. Wood <simon.wood@r-project.org>.
}

