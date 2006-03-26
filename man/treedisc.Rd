\name{treedisc}
\alias{treedisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prunes a level set tree or a tail tree
}
\description{
Prunes a level set tree or a tail tree so that it contains fewer levels.
The pruning makes plotting of the trees much faster.
The quality does not typically decrease essentially when a tree is pruned.
}
\usage{
treedisc(lst, pcf, ngrid = NULL, r = NULL, type = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{level set tree or a shape tree; output of for example "leafsfirst"}
  \item{pcf}{piecewise constant function; 
  output of for example "pcf.kern" or "pcf.func";
  "lst" should be the level set tree or shape tree of "pcf"}
  \item{ngrid}{positive integer; the number of levels which the pruned tree
  will have}
  \item{r}{vector of positive reals in increasing order; 
  we may give the set of levels explicitely with the argument "r";
  then the levels need not be equispaced as when we apply argument "ngrid"}
  \item{type}{"lst" or "shape"; not needed}
}
%\details{}
\value{
a level set tree or a shape tree which contains fewer levels than the input
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{leafsfirst}}
}
\examples{
dendat<-sim.data(n=100,type="mulmod")
pcf<-pcf.kern(dendat,h=1,N=c(32,32))
lst<-leafsfirst(pcf)          # level set tree
td<-treedisc(lst,pcf,ngrid=30)

plotvolu(td)

r<-c(0.01,0.03,0.04)
td<-treedisc(lst,pcf,r=r)

plotvolu(td)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS


