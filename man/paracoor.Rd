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
paracoor(X, Y=NULL, xmargin = 0.1,
paletti = matrix("black",dim(X)[1],1),
noadd = TRUE, verti = NULL, cex.axis = 1,
points=TRUE, col.verti="black", col.verti.y="red", digits=3,
arg=NULL, colarg="red", lwd=1, cex=1)
}
\arguments{
  \item{X}{n*d data matrix of real numbers}
  \item{Y}{n*1 data matrix of real numbers; response variable}
  \item{xmargin}{positive real number; empty space in left and right}
  \item{paletti}{n vector of colors; the lines will be colored with 
   these colors}
  \item{noadd}{TRUE or FALSE; if FALSE, then the lines will be added to an 
   existing plot}
  \item{verti}{real number or NULL; gives the x coordinate of a vertical line;
   can be used to divide the plot vertically to two parts} 
  \item{cex.axis}{the magnification to be used for axis annotation}
  \item{points}{TRUE or FALSE; whether points are plotted and not only lines}
  \item{col.verti}{a caharacter; the name of a color}
  \item{col.verti.y}{a caharacter; the name of a color}
  \item{digits}{a natural number; the number of digits in the annotation}
  \item{arg}{d*1 vector; the current value of the explanatory variable}
  \item{colarg}{a character; the name of a color}
  \item{lwd}{a positive real number; the width of the lines}
  \item{cex}{a positive real number; the largness of the circles}
}
\value{
Makes a plot on the graphics window
}
\references{Inselberg (1985), Wegman (1990)}
\author{Jussi Klemela}

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
X<-sim.data(n=100,type="mulmod")
paracoor(X)

set.seed(1)
X<-matrix(rnorm(300),100,3)
paracoor(X)

set.seed(1)
X<-matrix(rnorm(400),100,4)
paracoor(X)

}
\keyword{multivariate}


