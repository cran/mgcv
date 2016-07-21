\name{bug.reports.mgcv}
\alias{bug.reports.mgcv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Reporting mgcv bugs.}
\description{\code{mgcv} works largely because many people have reported bugs over the years. If you find something that looks like a bug, please report it, so that the package can be improved. \code{mgcv} does not have a large development budget, so it is a big help if bug reports follow the following guidelines.

The ideal report consists of an email to \email{simon.wood@r-project.org} with a subject line including \code{mgcv} somewhere, containing

\enumerate{
\item The results of running \code{\link{sessionInfo}} in the R session where the problem occurs. This provides platform details, R and package version numbers, etc.
\item A brief description of the problem.
\item Short cut and paste-able code that produces the problem, including the code for loading/generating the data (using standard R functions like \code{load}, \code{read.table} etc).
\item Any required data files. If you send real data it will only be used for the purposes of de-bugging.
}

Of course if you have dug deeper and have an idea of what is causing the problem, that is also helpful to know, as is any suggested code fix. (Don't send a fixed package .tar.gz file, however - I can't use this).   
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}


