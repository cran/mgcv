\name{uniquecombs}
\alias{uniquecombs}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{find the unique rows in a matrix }
\description{
This routine returns a matrix or data frame containing all the unique rows of the
matrix or data frame supplied as its argument. That is, all the duplicate rows are
stripped out. Note that the ordering of the rows on exit need not be the same
as on entry. It also returns an index attribute for relating the result back 
to the original matrix.
}
\usage{
uniquecombs(x,ordered=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{ is an \R matrix (numeric), or data frame. }
 \item{ordered}{ set to \code{TRUE} to have the rows of the returned object in the same order regardless of input ordering.}
}
\details{ Models with more parameters than unique combinations of
  covariates are not identifiable. This routine provides a means of
  evaluating the number of unique combinations of covariates in a
  model. 

  When \code{x} has only one column then the routine
  uses \code{\link{unique}} and \code{\link{match}} to get the index. When there are
  multiple columns then it uses \code{\link{paste0}} to produce labels for each row, 
  which should be unique if the row is unique. Then \code{unique} and \code{match} 
  can be used as in the single column case. Obviously the pasting is inefficient, but 
  still quicker for large n than the C based code that used to be called by this routine, which 
  had O(nlog(n)) cost. In principle a hash table based solution in C 
  would be only O(n) and much quicker in the multicolumn case. 
   
 \code{\link{unique}} and \code{\link{duplicated}}, can be used 
 in place of this, if the full index is not needed. Relative performance is variable. 

 If \code{x} is not a matrix or data frame on entry then an attmept is made to coerce 
 it to a data frame. 
}
\value{
A matrix or data frame consisting of the unique rows of \code{x} (in arbitrary order).

The matrix or data frame has an \code{"index"} attribute. \code{index[i]} gives the row of the returned 
matrix that contains row i of the original matrix. 
 
}

\seealso{\code{\link{unique}}, \code{\link{duplicated}}, \code{\link{match}}.}

\author{ Simon N. Wood \email{simon.wood@r-project.org} with thanks to Jonathan Rougier}


\examples{
require(mgcv)

## matrix example...
X <- matrix(c(1,2,3,1,2,3,4,5,6,1,3,2,4,5,6,1,1,1),6,3,byrow=TRUE)
print(X)
Xu <- uniquecombs(X);Xu
ind <- attr(Xu,"index")
## find the value for row 3 of the original from Xu
Xu[ind[3],];X[3,]

## same with fixed output ordering
Xu <- uniquecombs(X,TRUE);Xu
ind <- attr(Xu,"index")
## find the value for row 3 of the original from Xu
Xu[ind[3],];X[3,]


## data frame example...
df <- data.frame(f=factor(c("er",3,"b","er",3,3,1,2,"b")),
      x=c(.5,1,1.4,.5,1,.6,4,3,1.7))
uniquecombs(df)
}
\keyword{models} \keyword{regression}%-- one or more ..


