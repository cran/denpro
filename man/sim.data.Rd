\name{sim.data}
\alias{sim.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generates data for illustrative purposes}
\description{
Returns a random sample from some distributions,
to illustrate some visulization tools.
Returns also the density (as a piecewise constant function) 
for some examples, or the distribution function.
}
\usage{
sim.data(n = NULL, seed = 1, N = NULL, type = "mulmod",
M=NULL, sig=NULL, p=NULL, d=NULL,
cova=NULL, marginal="student", t=NULL, df=NULL, distr=FALSE)
}

\arguments{
  \item{n}{positive integer; size of the sample to be generated}
  \item{seed}{real number; seed for the random number generator.}
  \item{N}{2*1 vector of positive integers; the size of the grid where
   the piecewise constant function is evaluated }
  \item{type}{"mixt", "mulmod", "fox", "tetra3d", "penta4d", 
   "cross", "gauss", "student", "gumbel", "1d2modal", or "claw".}
  \item{M}{mixnum*d-matrix; rows of M are means of the Gaussians in
   the mixture. We have a mixture of "mixnum" Gaussians, whose
   dimension is d.}
  \item{sig}{mixnum*d-matrix; rows of sig are the diagonals of the
   covariance matrices of the mixtures.}
  \item{p}{mixnum-vector; weights for the members of the mixture.
   The sum of elements of "p" is 1.}
  \item{d}{positive integer; dimension of the vectors of the sample 
   to be generated, need to be given only when type="mixt" and d=1}
  \item{cova}{Covariance matrix for the Gauss or Student copulas}
  \item{marginal}{Marginal for the Gauss or Student copulas} 
  \item{t}{If marginal="student", gives the degrees of freedom}
  \item{df}{Degrees of freedom for the Student copula}
  \item{distr}{internal (implemented for "1d2modal")
   TRUE, if one wants the distribution function instead of the
   density function}
}
\details{
When type="mixt", generates data from a mixture of Gaussians.
When type="mulmod", the density is 3-modal.
When type="fox", the density has multimodal level sets.
}
\value{
If "n" is not NULL, then the function returns a n*d-data matrix or a
n*2-data matrix,
if "N" is not NULL, then the function returns a piecewise constant 
function on the grid of size N[1]*N[2],
if the both are NULL, then the function returns the mean, 
covariance, and the weights of the mixture components

}
\author{ Jussi Klemelä }

%\seealso{\code{\link{simmix}}}

\examples{
d<-2
mixnum<-3
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)   
M[2,]<-c(4,0)  
M[3,]<-c(0,4)   
sig<-matrix(1,mixnum,d)
p0<-1/mixnum
p<-p0*rep(1,mixnum)
n<-100
dendat<-sim.data(type="mixt",n=n,M=M,sig=sig,p=p,seed=1)    
plot(dendat)

dendat<-sim.data(n=100)
plot(dendat)

N<-c(20,20)
pcf<-sim.data(N=N)
dp<-draw.pcf(pcf,pnum=c(30,30))
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

sim.data()

type="fox"
dendat<-sim.data(n=100,type=type)
plot(dendat)

pcf<-sim.data(N=N,type=type)
dp<-draw.pcf(pcf,pnum=c(30,30))
contour(dp$x,dp$y,dp$z,drawlabels=FALSE)

}
\keyword{datagen}% at least one, from doc/KEYWORDS


