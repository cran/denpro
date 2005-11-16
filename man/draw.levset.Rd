\name{draw.levset}
\alias{draw.levset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a level set of a 2D function 
}
\description{
Plots a level set of a piecewise constant 2D function.
}
\usage{
draw.levset(pcf, lev, bary = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of "pcf.kern" or "eval.func"}
  \item{lev}{real number; gives the level of the level set}
  \item{bary}{vector of 2 reals; the barycenter; if given will be plotted}
}
%\details{}
\value{
Makes a plot to the graphics window
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemelä }
%\note{ ~~further notes~~ }

\seealso{
\code{\link{pcf.kern}},
\code{\link{eval.func}}
}
\examples{
dendat<-sim.exa(n=100)
pcf<-pcf.kern(dendat,h=1,N=c(25,25))

draw.levset(pcf,lev=0.01)

}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
