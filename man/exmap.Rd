\name{exmap}
\alias{exmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes a scale of excess mass profiles
}
\description{
Calculates a scale of excess mass profiles
from a scale of estimates given as level set trees
}
\usage{
exmap(estiseq, mt, hind = NULL, hseq = NULL, n = NULL, moteslist = NULL, 
ylist = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estiseq}{a scale of estimates, output for example of "lstseq.kern"}
  \item{mt}{mode graph; output of "modegraph"}
  \item{hind}{the indexes of the estimates for which the excess mass profiles
  will be calculated }
  \item{hseq}{ the sequence of smoothing parameters corresponding 
   to estimates }
  \item{n}{ sample size}
  \item{moteslist}{ internal }
  \item{ylist}{ internal, alternative to hind }
}
%\details{}
\value{
A list of excess mass profiles
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemela }
%\note{ ~~further notes~~ }

\seealso{
\code{\link{plotexmap}}, 
\code{\link{modegraph}}
}
\examples{
dendat<-sim.data(n=200,type="mulmod")

h1<-0.9
h2<-2.2
lkm<-5
hseq<-hgrid(h1,h2,lkm)

N<-c(16,16)
estiseq<-lstseq.kern(dendat,hseq,N,lstree=TRUE)

mt<-modegraph(estiseq)

horilines<-c(1,2,3,4)

sp<-exmap(estiseq,mt,horilines)

plotexmap(sp,mt,ylim=c(1,4.5))

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ smooth }% __ONLY ONE__ keyword per line
