\name{shape2d}
\alias{shape2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns a 2D volume function or 2D probability content function }
\description{
Calculates a 2D volume function or 2D probability content function
from a sequence of shape trees, corresponding to a sequence of level sets of a
multivariate function to be visualized.
A 2D volume function and 2D probability content function
are visualizations of a multivariate density with a 2D function.
One joins a series of radius functions or tail probability functions
to get a 2D visualization.  
Radius functions and tail probability functions are 1D functions
which visualize the shape of level sets of a multivariate density.
}
\usage{
shape2d(shtseq, gnum = 500, type = "radius", type2 = "slice", 
gnum2 = 1000, ngrid = 30, norma = FALSE, xmax = 10, modelim = 2, 
exmalim = NULL, maxnum = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shtseq}{A list of shape trees, made by function "stseq";
   the shape trees in the list correspond to a grid of level sets of the 
   function to be visualized}
  \item{gnum}{ Number of grid points in the radius function
   or in the tail probability function
   (used when one transforms a shape tree to a 1D function) }
  \item{type}{"radius" or "proba"; whether a 2D volume function or 2D
   probability content function will be calculated }
  \item{type2}{"slice" or "boundary"; whether the 1D functions are slices or
   level sets of the 2D function. The option "boundary" is not recommended
   at the moment }
  \item{gnum2}{when type2="boundary", then gnum*gnum is the size of the grid
   of the 2D function, and "gnum2" is the grid used for transforming
   shape trees to 1D functions }
  \item{ngrid}{ positive integer; one prunes the 1D functions to have "ngrid"
   level sets }
  \item{norma}{ TRUE if one uses dimension normalization for the volumes;
   then volumes are tranformed to (volume/V)pow(1/d), where V is the
   volume of the d-dimensional unit ball }
  \item{xmax}{ internal }
  \item{modelim}{ internal }
  \item{exmalim}{ internal }
  \item{maxnum}{ internal }
}
%\details{}
\value{
  A list containing elements x, y, z, whose meaning is the same as in 
  function "persp".
}
\references{ 
Jussi Klemela (2005).
Visualization of the spread of multivariate distributions.
}
\author{Jussi Klemela }
%\note{ ~~further notes~~ }
\seealso{
\code{\link{stseq}}, 
\code{\link{plotvolu2d}}
}
\examples{
N<-c(30,30)  
lnum<-20
func<-"gumbel"
marginal<-"normal"
g<-2
ver<-3
support<-c(-ver,ver,-ver,ver)
st<-stseq(N,lnum,func=func,marginal=marginal,g=g,support=support)

gnum<-50 
ngrid=50 
vd<-shape2d(st,gnum=gnum,ngrid=ngrid) 

plotvolu2d(vd)

type<-"proba"
vd3<-shape2d(st,gnum=gnum,type=type,ngrid=ngrid,norma=TRUE)

plotvolu2d(vd3)

}
\keyword{ multivariate }
