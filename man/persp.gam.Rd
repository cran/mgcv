\name{persp.gam}
\title{Perspective Plot of GAM objects}
\usage{
persp.gam(x, view=NULL, slice=list(), sizes=c(20,20), mask=FALSE, se=2,  
      theta=0, phi=15, r = sqrt(3), d = 1, scale = TRUE, expand = 1, col = NULL,
      border = NULL, ltheta = -135, lphi = 45, shade = 0.75, box = TRUE, axes = TRUE, 
      nticks = 5, ticktype = "detailed",...)   
}
\alias{persp.gam}
\arguments{
  \item{x}{a \code{gam} object, produced by \code{gam()}}
  \item{view}{a vector containing the names of the two covariates to be displayed on the 
   x and y dimensions of the plot.}
  \item{slice}{a named list of the values used for the covariates not directly displayed by 
   the diagram. Either in numbers or as functions of the two \code{view} parameters. Any 
   covariates not mentioned in the list will default to the mean of the gam's datapoints 
   in that dimension.}
  \item{sizes}{The number of grid squares in each direction used for calculating the 
   plotted surface.} 
 \item{mask}{If \code{FALSE} (default) the whole surface is displayed, otherwise those parts 
   of the surface far from all the datapoints are removed. Numeric values of \code{mask} 
   determine what is considered far. If mask is TRUE the maximum nearest neighbour distance
   for the projection of the datapoints onto each of the \code{view} dimensions is used.}
 \item{se}{if positive meshs this many standard errors above (red) and below (green) the 
  main (grey) function surface are added to the plot. Zero, negative or \code{FALSE} values suppress production of these extra surfaces}
\item{theta}{angle defining the azimuthal viewing direction.}
\item{phi}{the colatitudeof the viewing direction.}
\item{r}{the distance of the eyepoint from the centre of the plotting box.}
\item{d}{a value which can be used to vary the strength of
    the perspective transformation.  Values of \code{d} greater
    than 1 will lessen the perspective effect and values less
    and 1 will exaggerate it.}
  \item{scale}{before viewing the x, y and z coordinates of the
    points defining the surface are transformed to the interval
    [0,1].  If \code{scale} is \code{TRUE} the x, y and z coordinates
    are transformed separately.  If \code{scale} is \code{FALSE}
    the coordinates are scaled so that aspect ratios are retained.
    This is useful for rendering things like DEM information.}
  \item{expand}{a expansion factor applied to the \code{z}
    coordinates. Often used with \code{0 < expand < 1} to shrink the
    plotting box in the \code{z} direction.}
 \item{col}{the color of the surface facets.}
  \item{border}{the color of the line drawn around the surface facets.
    A value of \code{NA} will disable the drawing of borders.  This is
    sometimes useful when the surface is shaded.}
  \item{ltheta}{if finite values are specified for \code{ltheta}
    and \code{lphi}, the surface is shaded as though it was being
    illuminated from the direction specified by azimuth \code{ltheta}
    and colatitude \code{lphi}.}
  \item{lphi}{see ltheta.}
  \item{shade}{the shade at a surface facet is computed as
    \code{((1+d)/2)^shade}, where \code{d} is the dot product of
    a unit vector normal to the facet and a unit vector in the
    direction of a light source.  Values of \code{shade} close
    to one yield shading similar to a point light source model
    and values close to zero produce no shading.  Values in the
    range 0.5 to 0.75 provide an approximation to daylight
    illumination.}
  \item{box}{should the bounding box for the surface be displayed.
    The default is \code{TRUE}.}
  \item{axes}{should ticks and labels be added to the box.  The 
    default is \code{TRUE}.  If \code{box} is \code{FALSE} then no
    ticks or labels are drawn.}
 \item{nticks}{the (approximate) number of tick marks to draw on the
    axes.  Has no effect if \code{ticktype} is "simple".}
  \item{ticktype}{character: "simple" draws just an arrow parallel to 
    the axis to indicate direction of increase; "detailed" draws normal
    ticks as per 2D plots.}
  \item{...}{ other arguments.}
}
\value{ If \code{se=FALSE} the matrix of predicted values for the gam object is returned.
   This can then be rotated and manipulated in \code{persp} without needing to 
   reconstruct the surface each time.Otherwise a three element list containing the 
   predicted values and the upper and lower confidence intervals is returned.}

\description{
  This function draws a perspective plot of a slice through the response 
  surface of a gam object. Meshes indicating confidence intervals around the 
  surface can also be displayed and most of the parameters in \code{persp} 
  set within it.
}
\details{
  A grid of values covering the projection of the range of data on which the \code{gam} 
  object, \code{x}, is based are defined and the predicted values for each value them found. 
  A perspective plot of this is then drawn. The standard errors around the surface 
  can be indicated, and parts of the surface far from the data removed, if desired.
}

\author{ Mike Lonergan \email{mel@st-and.ac.uk}}

\seealso{
  \code{\link{persp}} and \code{\link{gam}}.
}
\examples{
library(mgcv)
set.seed(c(1,2,3))
n<-200
sig2<-4
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
y<-x0^2+x1*x2 +runif(n,-0.3,0.3)
g<-gam(y~s(x0,x1,x2))
persp.gam(g,se=FALSE)  # displays the prediction surface in x0, x1 
     # with x2 fixed at its mean and no confidence intervals
persp.gam(g, view=c("x1","x2"),slice=list(x0=0.75)) # a different slice
persp.gam(g, view=c("x0","x2"),slice=list(x1="x2"))->p 
                       # a view of another slice is drawn and the surfaces saved 
persp(p[[1]],theta=60)          # another view of the main surface above
persp.gam(g, view=c("x1","x2"),slice=list(x0=0.75),theta=210,phi=40,mask=0.07)
                       # areas where there is no data not plotted
}
\keyword{models} \keyword{smooth} \keyword{regression}