\name{branchmap}
\alias{branchmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes a branching map from a sequence of level set trees}
\description{
Branching map visualizes the levels of branching of level set trees
of estimates belonging to a scale of estimates.
It visualizes also the excess masses of the roots of the branches. 
}
\usage{
branchmap(estiseq, hseq = NULL, levnum = 80, paletti = NULL, 
rootpaletti = NULL, type = "jump")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estiseq}{A sequence of estimates and level set trees of the estimates.
   Output of function lstseq.kern or function lstseq.carthisto.}
  \item{hseq}{The sequence of smoothing parameters of the scale of estimates.}
  \item{levnum}{The number of level sets used to approximate the 
   level set trees.}
  \item{paletti}{A sequence of color names; colors for each branch,
   other than the root branches.}
  \item{rootpaletti}{A sequence of color names; colors for the root branches.}
  \item{type}{internal}
}
%\details{}
\value{
  A representation as a list of a 2D function
  \item{level}{x-coordinate is the level of the level sets}
  \item{h}{y-coordinate is the smoothing parameter}
  \item{z}{z-coordinate is the excess mass}
  \item{col}{colors for the graph of the 2D function}
}
%\references{http://denstruct.net}
\author{Jussi Klemela}
%\note{ ~~further notes~~ }

\seealso{ 
\code{\link{lstseq.kern}},
\code{\link{plotbranchmap}}
}
\examples{
dendat<-sim.data(n=200,type="mulmod")

h1<-0.9
h2<-2.2
lkm<-5
hseq<-hgrid(h1,h2,lkm)

N<-c(16,16)
estiseq<-lstseq.kern(dendat,hseq,N,lstree=TRUE)

bm<-branchmap(estiseq)

plotbranchmap(bm)

}

\keyword{smooth}  
\keyword{multivariate}  % at least one, from doc/KEYWORDS



