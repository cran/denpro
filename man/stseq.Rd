\name{stseq}
\alias{stseq}
\title{ Calculates a sequence of radius functions
from a sequence of level sets }
\description{
 Calculates a sequence of radius functions from level sets of a function
 of a given parametric form.
}
\usage{
stseq(N, lnum, refe = NULL, func = NULL, dendat = NULL, h = NULL, Q = NULL, 
kernel = "epane", hw = NULL, sig = rep(1, length(N)), support = NULL, 
theta = NULL, M = NULL, p = NULL, mul = 3, t = rep(1, length(N)), 
marginal = "normal", r = 0, mu = NULL, xi = NULL, Omega = NULL, alpha = NULL, 
df = NULL, g = 1, base = 10)
}
\arguments{
  \item{N}{d-vector of integers >1; the size of the grid where the function is
  evaluated}
  \item{lnum}{positive integer; the number of level sets from which the
  transforms are calculated  }
  \item{refe}{d-vector of real numbers; the reference point for the shape
   trees; if refe=NULL, then the location of the maximum of the function will
   be used}
  \item{func}{character; the name of the function, see the "pcf.func" for the
  possibilities}
  \item{dendat}{n*d matrix of data; when the function is a kernel estimate one
  gives the data as argument}
  \item{h}{postive real number; the smoothing parameter of the kernel estimate}
  \item{Q}{ internal }
  \item{kernel}{"gauss" or "epane"; the kernel of the kernel estimate }
  \item{hw}{ positive real number; the time window of the kernel estimate,
  used only for time series data }
  \item{sig}{ see "pcf.func" }
  \item{support}{ see "pcf.func" }
  \item{theta}{ see "pcf.func" }
  \item{M}{ see "pcf.func" }
  \item{p}{ see "pcf.func" }
  \item{mul}{ see "pcf.func" }
  \item{t}{ see "pcf.func" }
  \item{marginal}{ see "pcf.func" }
  \item{r}{ see "pcf.func" }
  \item{mu}{ see "pcf.func" }
  \item{xi}{ see "pcf.func" }
  \item{Omega}{ see "pcf.func" }
  \item{alpha}{ see "pcf.func" }
  \item{df}{ see "pcf.func" }
  \item{g}{ see "pcf.func" }
  \item{base}{ positve integer or NULL; the base of the logarithm, when the
  logarithmic spacing is used for the level sets, if "base" is NULL, then the
  level sets are equispaced}
}
%\details{}
\value{
 A list of radius transforms.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }

\seealso{ 
\code{\link{shape2d}}, 
\code{\link{plotvolu2d}} 
}
\examples{
N<-c(32,32)
lnum<-30
func<-"prod"
marginal<-"student"
yla<-5
support<-c(-yla,yla,-yla,yla)
g<-c(1,3)
st<-stseq(N,lnum,func=func,marginal=marginal,support=support,g=g)

gnum<-70  
vd<-shape2d(st,gnum=gnum)
plotvolu2d(vd)

}
\keyword{multivariate}



