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
paracoor(dendat)
}
\arguments{
  \item{dendat}{n*d data matrix of real numbers}
}
\value{
Makes a plot on the graphics window
}
\references{Inselberg (1985), Wegman (1990)}
\author{Jussi Klemelä}

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
dendat<-sim.exa(n=100)
paracoor(dendat)

set.seed(1)
dendat<-matrix(rnorm(300),100,3)
paracoor(dendat)

set.seed(1)
dendat<-matrix(rnorm(400),100,4)
paracoor(dendat)

}
\keyword{multivariate}
\keyword{hplot}
