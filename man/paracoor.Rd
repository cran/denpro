\name{paracoor}
\alias{paracoor}
\title{Makes a parallel coordinate plot}
\description{
Makes a parallel coordinate plot of a data matrix:
the coordinates of the observations are plotted so that
the coordinate axes are drawn parallel to each other and the
coordinates of a single observation are joined with lines.
}
\usage{
paracoor(dendat, xmargin = 0.1,
paletti = matrix("black",dim(dendat)[1],1),
noadd = TRUE, verti = NULL, cex.axis = 1)
}
\arguments{
  \item{dendat}{n*d data matrix of real numbers}
  \item{xmargin}{positive real number; empty space in left and right}
  \item{paletti}{n vector of colors; the lines will be colored with 
   these colors}
  \item{noadd}{TRUE or FALSE; if FALSE, then the lines will be added to an 
   existing plot}
  \item{verti}{real number or NULL; gives the x coordinate of a vertical line;
   can be used to divide the plot vertically to two parts} 
  \item{cex.axis}{the magnification to be used for axis annotation}
}
\value{
Makes a plot on the graphics window
}
\references{Inselberg (1985), Wegman (1990)}
\author{Jussi Klemela}

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
dendat<-sim.data(n=100,type="mulmod")
paracoor(dendat)

set.seed(1)
dendat<-matrix(rnorm(300),100,3)
paracoor(dendat)

set.seed(1)
dendat<-matrix(rnorm(400),100,4)
paracoor(dendat)

}
\keyword{multivariate}

