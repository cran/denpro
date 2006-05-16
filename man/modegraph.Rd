\name{modegraph}
\alias{modegraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Calculates a mode graph from a scale of estimates 
}
\description{
Calculates a mode graph from a scale of estimates.
The estimates in the scale are given as level set trees.
}

\usage{
modegraph(estiseq, hseq = NULL, paletti = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estiseq}{a scale of estimates, output for example of "lstseq.kern" }
  \item{hseq}{ the corresponding scale of smoothing parameters }
  \item{paletti}{ vector of color names }
}
%\details{}
\value{
mode graph
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{lstseq.kern}}, 
\code{\link{plotmodet}}
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

plotmodet(mt)

}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{ smooth }% __ONLY ONE__ keyword per line
