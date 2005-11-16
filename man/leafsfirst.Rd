\name{leafsfirst}
\alias{leafsfirst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates a level set tree, shape tree, or tail tree
}
\description{
Calculates a level set tree or shape tree from a piecewise constant
function, or calculates a tail tree from a data set.
}
\usage{
leafsfirst(pcf = NULL, lev = NULL, refe = NULL, type = "lst", 
levmet = "radius", ordmet = "etaisrec", ngrid = NULL, dendat = NULL, rho = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of for example "pcf.kern", or
  "eval.func"}
  \item{lev}{positive real number; a level of the level set whose shape tree
  will be calculated}
  \item{refe}{vector of d real numbers; the reference point of the shape tree}
  \item{type}{"lst" or "shape"; not needed}
  \item{levmet}{"radius" or "proba"; radius plot or probability content plot}
  \item{ordmet}{"etaisrec" or "baryrec"; a distance to a rectangle is defined
  to be the distance to the boundary or the distance to the barycenter;
  concerns shape trees}
  \item{ngrid}{a positive integer; the tree will be pruned to contain "ngrid" levels}
  \item{dendat}{n*d data matrix; gives the data set when a tail tree will be calculated}
  \item{rho}{positive real number; the resolution threshold for the tail tree}
}
%\details{}
\value{a level set tree, shape tree, or tail tree}
%\references{}
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{pcf.kern}},
\code{\link{treedisc}},
\code{\link{plotvolu}}, 
\code{\link{plotbary}}, 
}
\examples{
dendat<-sim.exa(n=100)
pcf<-pcf.kern(dendat,h=1,N=c(25,25))

# level set tree

lst<-leafsfirst(pcf)          

td<-treedisc(lst,pcf,ngrid=30)
plotvolu(td)

# shape tree

lev<-0.1*max(pcf$value)       # 10\% level set
st<-leafsfirst(pcf,lev=lev)   

tdst<-treedisc(st,pcf,ngrid=30)
plotvolu(tdst)

# tail tree

rho<-0.65
tt<-leafsfirst(dendat=dendat,rho=rho)

plotbary(tt,modelabel=FALSE)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS

