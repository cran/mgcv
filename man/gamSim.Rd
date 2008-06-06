\name{gamSim}
\alias{gamSim}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Simulate example data for GAMs}
\description{ Function used to simulate data sets to illustrate the use of 
\code{\link{gam}} and \code{\link{gamm}}. Mostly used in help files to keep down 
the length of the example code sections.

}

\usage{gamSim(eg=1,n=400,dist="normal",scale=2)}

%- maybe also `usage' for other objects documented here.
\arguments{
\item{eg}{ numeric value specifying the example required.}

\item{n}{ number of data to simulate.}

\item{dist}{character string which may be used to spcify the distribution of
the response.}

\item{scale}{Used to set noise level.}

  }

\details{See the source code for exactly what is simulated in each case. 
\item{1}{Gu and Wahba 4 univariate term example.}
\item{2}{A smooth function of 2 variables.}
\item{3}{Example with continuous by variable.}
\item{4}{Example with factor by variable.}
\item{5}{An additive example plus a factor variable.}
\item{6}{Additive + random effect.}
\item{7}{As 1 but with correlated covariates.}
}

\value{ Depends on \code{eg}, but usually a dataframe, which may also contain some information on 
the underlying truth. Sometimes a list with more items, including a data frame for model fitting.
See source code or helpfile examples where the function is used for further information.}



\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\seealso{ \code{\link{gam}}, \code{\link{gamm}}}

\examples{
## see ?gam
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









