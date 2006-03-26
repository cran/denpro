\name{pcf.kern}
\alias{pcf.kern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates a multivariate kernel estimate
}
\description{
Calculates a multivariate kernel estimate and gives the output as a piecewise
constant function object.
}
\usage{
pcf.kern(dendat, h, N, kernel = "epane", hw = NULL)
}
\arguments{
  \item{dendat}{n*d matrix of real numbers; the data matrix}
  \item{h}{positive real; smoothing parameter}
  \item{N}{vector of d positive dyadic integers; 
  the dimension of the grid where the
  kernel estimate will be evaluated; we evaluate the estimate on a regular
  grid which contains the support of the kernel estimate }
  \item{kernel}{"epane" or "gauss"; the kernel is either the
  Bartlett-Epanechnikov product kernel or the standard Gaussian}
  \item{hw}{positive integer; parameter for time localized kernel estimation;
  gives the smoothing prameter for the temporal smoothing}
}
%\details{}
\value{
a piecewise constant function object, see the web site
}
\references{http://www.rni.helsinki.fi/~jsk/denpro/}
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

 \seealso{
\code{\link{draw.pcf}}
\code{\link{lstseq.kern}}
}

\examples{
dendat<-sim.data(n=100,type="mulmod")

pcf<-pcf.kern(dendat,h=1,N=c(32,32))

dp<-draw.pcf(pcf)
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

pcf<-pcf.kern(dendat,h=1,N=c(32,32),hw=60)

dp<-draw.pcf(pcf)
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line

