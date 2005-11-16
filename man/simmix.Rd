\name{simmix}
\alias{simmix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generates a random sample from a mixture of Gaussians}
\description{
Generates a random sample from a mixture of multivariate Gaussians.
The Gaussians have diagonal covariance matrices.
}
\usage{
simmix(n, M, sig, p, seed, d=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{n}{positive integer; size of the sample to be generated}
 \item{M}{mixnum*d-matrix; rows of M are means of the Gaussians in
  the mixture. We have a mixture of "mixnum" Gaussians, whose
  dimension is d.
  }
 \item{sig}{mixnum*d-matrix; rows of sig are the diagonals of the
  covariance matrices of the mixtures.
  }
 \item{p}{mixnum-vector; weights for the members of the mixture.
  The sum of elements of "p" is 1.}
 \item{seed}{real number; seed for the random number generator.}
 \item{d}{positive integer; dimension of the vectors of the sample 
  to be generated, need to be given when d=1}
}
%\details{}
\value{n*d-matrix: n samples from a d-dimensional distribution}
\references{http://www.rni.helsinki.fi/~jsk/denpro}
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

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
n<-10
dendat<-simmix(n,M,sig,p,seed=1)    

}

\keyword{datagen}%-- one or more ...


