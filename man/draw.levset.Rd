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
draw.levset(pcf, lev=NULL, bary = NULL, propor=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of "pcf.kern" or "pcf.func"}
  \item{lev}{real number; gives the level of the level set}
  \item{bary}{vector of 2 reals; the barycenter; if given will be plotted}
  \item{propor}{0<propor<1; the level set whose level is "propor"
   times the maximum value of the function will be drawn}
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
\code{\link{pcf.func}}
}
\examples{
dendat<-sim.data(n=100,type="mulmod")
pcf<-pcf.kern(dendat,h=1,N=c(32,32))

draw.levset(pcf,lev=1.331979e-02)

draw.levset(pcf,propor=0.4)

}
\keyword{hplot}% __ONLY ONE__ keyword per line
