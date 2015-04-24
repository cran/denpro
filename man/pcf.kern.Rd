\name{pcf.kern}
\alias{pcf.kern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes a multivariate kernel estimate
}
\description{
Computes a multivariate kernel estimate and gives the output as a piecewise
constant function object.
}
\usage{
pcf.kern(dendat, h, N, kernel = "gauss", weights = NULL, support = NULL,
lowest = 0, radi = 0)
}
\arguments{
  \item{dendat}{n*d matrix of real numbers; the data matrix}
  \item{h}{d vector of positive real numbers; vector of smoothing parameters}
  \item{N}{vector of d positive dyadic integers; 
   the dimension of the grid where the
   kernel estimate will be evaluated; we evaluate the estimate on a regular
   grid which contains the support of the kernel estimate}
  \item{kernel}{"gauss", "epane", "bart", or "uniform"; 
   the kernel is either the standard Gaussian,
   Epanechnikov product kernel, Bartlett kernel, or uniform kernel}
  \item{weights}{n vector of nonnegative weights, where n is the sample size;
   sum of the elements of "weights" should be one;
   these are the weights of a time localized kernel estimator}
  \item{support}{2*d vector of reals gives the d intervals of a 
   rectangular support in the form c(low1,upp1,...,lowd,uppd)}
  \item{lowest}{a real value; the density estimate will be truncated to take
   value zero, if the value of the estimate is less or equal to "lowest"}
  \item{radi}{a nonnegative real number; the support is estimated as the smallest 
   rectangle containing the observations with an additional band whose
   width is equal to "radi"}
}
%\details{}
\value{
a piecewise constant function object, see the web site
}
\references{http://www.rni.helsinki.fi/~jsk/denpro/}
\author{Jussi Klemela}
%\note{ ~~further notes~~ }

 \seealso{
\code{\link{draw.pcf}}
\code{\link{lstseq.kern}}
}

\examples{
n<-100
dendat<-sim.data(n=n,type="mulmod")

h<-1
pcf<-pcf.kern(dendat,h=h,N=c(32,32))
dp<-draw.pcf(pcf)
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

d<-2
h<-(4/(d+2))^(1/(d+4))*n^(-1/(d+4))*apply(dendat,2,sd)
pcf<-pcf.kern(dendat,h=h,N=c(32,32))
dp<-draw.pcf(pcf)
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

# we use now nonuniform weighting of kernels

weights<-matrix(0,n,1)
threshold<-4
for (i in 1:n){
    eta<-(n-i)
    if (eta/h>threshold) result<-0 else result<-exp(-eta^2/(2*h^2))
    weights[i]<-result
}
weights<-weights/sum(weights)

pcf<-pcf.kern(dendat,h=1,N=c(32,32),weights=weights)

dp<-draw.pcf(pcf)
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line

