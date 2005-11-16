\name{sim.exa}
\alias{sim.exa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generates illustrative examples}
\description{
Returns a density (as a piecewise constant function) or a
random sample from the density,
for some illustrative purposes.
}
\usage{
sim.exa(n = NULL, seed = 1, N = NULL, type = "mulmod")
}

\arguments{
  \item{n}{sample size}
  \item{seed}{random seed}
  \item{N}{2*1 vector of positive integers; the size of the grid where
  the piecewise constant function is evaluated }
  \item{type}{"mulmod" or "fox"}
}
\details{
When type=="mulmod", the density is 3-modal,
when type="fox", the density has multimodal level sets
}
\value{
If "n" is not NULL, then the function returns a n*2-data matrix,
if "N" is not NULL, then the function returns a piecewise constant 
function on the grid of size N[1]*N[2],
if the both are NULL, then the function returns the mean, 
covariance, and the weights of the mixture components

}
\author{ Jussi Klemelä }

\seealso{\code{\link{simmix}}}
\examples{
dendat<-sim.exa(n=100)
plot(dendat)

N<-c(20,20)
pcf<-sim.exa(N=N)
dp<-draw.pcf(pcf,pnum=c(30,30))
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

sim.exa()

type="fox"
dendat<-sim.exa(n=100,type=type)
plot(dendat)

pcf<-sim.exa(N=N,type=type)
dp<-draw.pcf(pcf,pnum=c(30,30))
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

}
\keyword{datagen}% at least one, from doc/KEYWORDS



