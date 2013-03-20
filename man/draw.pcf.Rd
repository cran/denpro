\name{draw.pcf}
\alias{draw.pcf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepares the plotting of a 2D or 1D piecewise constant function
}
\description{
Calculates x, y, and z arguments for a piecewise constant function,
to be used in "persp" or "contour" functions,
or the "x" and "y" arguments to be used in "plot".
}
\usage{
draw.pcf(pcf, pnum = rep(32,length(pcf$N)), corona = 5, dens = FALSE, 
minval = 0, drawkern=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of "pcf.kern" or "pcf.func"}
  \item{pnum}{vector of 2 positive integers; dimension of the grid where the
   function will be plotted; not needed in the 1D case}
  \item{corona}{positive integer; gives the number of zeros around the support
   of the function; the plots look better when there are a corona of zeros
   around the function}
  \item{dens}{TRUE for densities, used in the 1D case instead of the corona}
  \item{minval}{real number, gives the level of the basis plane of the picture,
   one takes "minval" typically to be the minimum value of the function
   (for densities it would be zero)}
  \item{drawkern}{TRUE or FALSE; internal parameter; when one wants to use 
   function "drawkern" for calculation, then one should specify drawkern=TRUE,
   this can be used in the case that "pcf" is a kernel density estimate}
}
%\details{}
\value{
in the 2D case list of vectors "x" and "y" and matrix "z",
in the 1D case list of vectors "x" and "y"
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemela}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{pcf.kern}}, 
\code{\link{pcf.func}},
\code{\link{persp}},
\code{\link{contour}},   
}

\examples{
# 2D case
dendat<-sim.data(n=100,type="mulmod")
pcf<-pcf.kern(dendat,h=1,N=c(32,32))

dp<-draw.pcf(pcf,pnum=c(30,30))

contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

persp(dp$x,dp$y,dp$z)

# 1D case
set.seed(1)
dendat<-matrix(rnorm(20),20)   #10*1 data-matrix
pcf<-pcf.kern(dendat,h=1,N=c(25))

dp<-draw.pcf(pcf)

plot(dp$x,dp$y,type="l")

}
\keyword{hplot}% __ONLY ONE__ keyword per line

