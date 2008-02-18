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
levmet = "radius", ordmet = "etaisrec", ngrid = NULL, dendat = NULL, rho = 0,
propor= NULL, lowest = "dens")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of for example "pcf.kern", or
  "pcf.func"}
  \item{lev}{positive real number; a level of the level set whose shape tree
   will be calculated}
  \item{refe}{vector of d real numbers; the reference point of the shape tree}
  \item{type}{"lst" or "shape"; not needed}
  \item{levmet}{"radius" or "proba"; radius plot or probability content plot}
  \item{ordmet}{"etaisrec" or "baryrec"; a distance to a rectangle is defined
  to be the distance to the boundary or the distance to the barycenter;
  concerns shape trees}
  \item{ngrid}{a positive integer; the tree will be pruned to contain 
   "ngrid" levels}
  \item{dendat}{n*d data matrix; gives the data set when a tail tree 
   will be calculated}
  \item{rho}{positive real number; the resolution threshold for the tail tree}
  \item{propor}{0<propor<1; a shape tree of the level set whose level 
   is "propor" times the maximum value of the function will be calculated}
  \item{lowest}{a character; if lowest="dens", then it is assumed that
   the function whose level set tree we calculate is nonnegative}
}
%\details{}
\value{a level set tree, shape tree, or tail tree}
%\references{}
\author{Jussi Klemela}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{pcf.kern}},
\code{\link{treedisc}},
\code{\link{plotvolu}}, 
\code{\link{plotbary}}, 
}
\examples{
dendat<-sim.data(n=100,type="mulmod")
pcf<-pcf.kern(dendat,h=1,N=c(32,32))

# level set tree

lst<-leafsfirst(pcf)          

td<-treedisc(lst,pcf,ngrid=30)
plotvolu(td)

# shape tree

st<-leafsfirst(pcf,propo=0.1)   # 10\% level set  

tdst<-treedisc(st,pcf,ngrid=30)
plotvolu(tdst)

# tail tree

rho<-0.65
tt<-leafsfirst(dendat=dendat,rho=rho)

plotbary(tt)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{smooth}

