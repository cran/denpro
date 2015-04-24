\name{findbnodes}
\alias{findbnodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finds pointers to the nodes of a level set tree
}
\description{
Finds pointers to the nodes of a level set tree.
The nodes are such that their excess mass is among "modenum" largest excess masses, 
and they start a new branch.
}
\usage{
findbnodes(lst, modenum=1, num=NULL)
}
\arguments{
  \item{lst}{level set tree; a list of vectors.
   The list contains at least vectors "level", "volume", and "parent". 
   For example, function "leafsfirst" computes a level set tree.}
  \item{modenum}{positive integer; the number of nodes whose pointers are computed}
  \item{num}{a technical parameter}
}
%\details{}
\value{
a vector of integers; pointers to the nodes of the level set tree
}
%\references{http://www.rni.helsinki.fi/~jsk/denpro/}
\author{Jussi Klemela}
%\note{ ~~further notes~~ }

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
N<-rep(20,d)
pcf<-pcf.kern(dendat,h,N)
lst<-leafsfirst(pcf)

nodes<-findbnodes(lst,modenum=3)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line

