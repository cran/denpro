\name{plotvolu2d}
\alias{plotvolu2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Makes a perspective plot of 
a 2D volume function or a 2D probability content function}
\description{
Plots the output of function "shape2d";
this is a 2D volume function or a 2D probability content function.
}
\usage{
plotvolu2d(vd, theta = NULL, phi = NULL, typer = "flat")
}
\arguments{
  \item{vd}{Output of function "shape2d";
   a list containing fields x, y, z; the meaning of the fealds is the same as
   in the function "persp" }
  \item{theta}{ Viewing angle (left-right); the same parameter as in "persp"  }
  \item{phi}{ Viewing angle (up-down), the same as in "persp" }
  \item{typer}{"flat" or "dens";  internal }
}
%\details{}
\value{Plots a perspective plot in the graphics window}
\references{ 
Jussi Klemela (2005).
Visualization of the spread of multivariate distributions.}
\author{Jussi Klemela}
%\note{ }
\seealso{\code{\link{shape2d}}}
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

plotvolu2d(vd,theta=30)

}
\keyword{ hplot }



