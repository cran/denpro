\name{tree.segme}
\alias{tree.segme}
\title{ Returns the segmentation of the nodes of a visualization tree }
\description{
Returns the segmentation of the nodes of a visualization tree.
When the visualization tree is a tail tree,
finds a grouping for the observations.
The grouping may be used to enhance scatter plots,
graphical matrices, and parallel coordinate plots.
}
\usage{
tree.segme(tt, paletti = NULL, pcf = NULL)
}
\arguments{
  \item{tt}{ visualization  tree;
   for example a tail tree }
  \item{paletti}{a sequence of numbers or colors (group labels)}
  \item{pcf}{internal; piecewise constant function}
}
%\details{}
\value{
Vector of length n, where n is the number of nodes in the visualization
tree.
The elements of the vector give the segmentation labels.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{paracoor}}
}
\examples{
dendat<-sim.data(n=1000,type="cross",seed=1)
rho<-1.1
tt<-leafsfirst(dendat=dendat,rho=rho)

ts<-tree.segme(tt)

plot(dendat,col=ts)

paracoor(dendat,paletti=ts)

}
\keyword{ cluster }% at least one, from doc/KEYWORDS

