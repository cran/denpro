\name{plotexmap}
\alias{plotexmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Plots a scale of excess mass profiles
}
\description{
Plots a scale of excess mass profiles.
}
\usage{
plotexmap(sp, mt, xaxt = "n", lift = 0.1, leaflift = 0.1, ylim = NULL, 
leafcolors = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp}{scale of excess mass profiles; output of "exmap" }
  \item{mt}{mode graph; output of "modegraph" }
  \item{xaxt}{"n" or "y"; whether x-axis will be plotted }
  \item{lift}{positive real number; space between the lines of the excess
  mass profiles }
  \item{leaflift}{ positive real number; lifting for the bullets associated
  with the leaf nodes of the excess mass profiles }
  \item{ylim}{ 2-vector; gives the interval of the y-axis}
  \item{leafcolors}{a paletti of the colors for the labeling of the 
   leaf nodes }
}
%\details{}
\value{
Plots the scale of excess mass profiles to the graphic device
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{ 
\code{\link{exmap}}, 
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
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{smooth }% __ONLY ONE__ keyword per line
\keyword{hplot }
