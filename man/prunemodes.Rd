\name{prunemodes}
\alias{prunemodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prunes modes away from a level set tree or a shape tree}
\description{
Prunes smallest modes in terms of the excess mass away, from a level set tree
or froma a shape tree.
The tree is interpreted as a 1D function (in the sense of a volume function,
radius function, or a tail probability function).
The pruning may be useful to get cleaner plots,
especially a location plot of a shape tree may suffer from spurious modes.
}
\usage{
prunemodes(lst, modenum = 1, num = NULL, exmalim = NULL, maxnum = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{A level set tree or a shape tree}
  \item{modenum}{positive integer; the desired number of modes}
  \item{num}{ internal }
  \item{exmalim}{ internal }
  \item{maxnum}{ internal }
}
%\details{}
\value{
A level set tree or a shape tree,
with an additional field "exma.of.modes".
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{leafsfirst}},
\code{\link{treedisc}},
}
\examples{
dendat<-sim.data(n=100,type="cross")
pcf<-pcf.kern(dendat,h=1,N=c(32,32))
st<-leafsfirst(pcf,propor=0.01)   

td<-treedisc(st,pcf,ngrid=60)
plotbary(td)

td2<-prunemodes(td,modenum=4)
plotbary(td2)


}
\keyword{hplot}


