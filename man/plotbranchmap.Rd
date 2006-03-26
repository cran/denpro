\name{plotbranchmap}
\alias{plotbranchmap}
\title{Plots a branching map}
\description{
Plots a branching map to the graphics window.
}
\usage{
plotbranchmap(bm, phi = 55, theta = 30)
}
\arguments{
  \item{bm}{A branching map; output of function "branchmap".}
  \item{phi}{The phi-parameter of function "persp".}
  \item{theta}{The theta-parameter of function "persp".}
}
\details{Function "plotbranchmap" just calls function "persp".}
\value{
Makes a plot to the current graphics window.
}
\references{http://denstruct.net}
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{branchmap}}
}
\examples{
dendat<-sim.data(n=200,type="mulmod")

h1<-0.9
h2<-2.2
lkm<-5
hseq<-hgrid(h1,h2,lkm)

N<-c(16,16)
estiseq<-lstseq.kern(dendat,hseq,N)

bm<-branchmap(estiseq)

plotbranchmap(bm)

## The function is currently defined as
function(bm,phi=55,theta=30)
{

persp(x=bm$level,y=bm$h,z=bm$z, 
xlab="level",ylab="h",zlab="excess mass",
ticktype="detailed",
col=bm$col,phi=phi,theta=theta) 

  }
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
