\name{leafsfirst.adagrid}
\alias{leafsfirst.adagrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes a level set tree of a discretized kernel estimator with an adaptive grid 
}
\description{
Computes a level set tree of a discretized kernel estimator with an adaptive grid .
}
\usage{
leafsfirst.adagrid(pcf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of "pcf.greedy.kernel"}
}
%\details{}
\value{a level set tree}
\references{http://jussiklemela.com/Art/volume}
\author{Jussi Klemela}
\note{Function pcf.greedy.kernel is in package "delt"}

\seealso{
\code{\link{leafsfirst}}
}
\examples{
# generate data
seed<-1
n<-50
d<-2
l<-3; D<-4; c<-D/sqrt(2)
M<-matrix(0,l,d); M[2,]<-c; M[3,]<--c
sig<-matrix(1,l,d)
p<-rep(1/l,l)
dendat<-sim.data(type="mixt",n=n,M=M,sig=sig,p=p,seed=seed)

h<-(4/(d+2))^(1/(d+4))*n^(-1/(d+4))*apply(dendat,2,sd)
minobs<-3
#pcf<-pcf.greedy.kernel(dendat,h,minobs=minobs,type="greedy")
#lst<-leafsfirst.adagrid(pcf)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{smooth}


