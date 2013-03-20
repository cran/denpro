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
draw.levset(pcf, lev=NULL, bary = NULL, propor = 0.1, col=NULL, 
bound = NULL, dendat = NULL, xaxt = "s", yaxt = "s", cex.axis = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of "pcf.kern" or "pcf.func"}
  \item{lev}{real number; gives the level of the level set}
  \item{bary}{vector of 2 reals; the barycenter; if given will be plotted}
  \item{propor}{0<propor<1; the level set whose level is "propor"
   times the maximum value of the function will be drawn}
  \item{col}{the color of the lines of the plot;
   for example "red"}
  \item{bound}{vector of 4 real numbers; c(xmin,xmax,ymin,ymax);
   gives the values for the parameters "xlim" and "ylim" in the "plot" 
   function; if "bound" is not specified, then the support of function "pcf"
   will be used}
  \item{dendat}{n*2-matrix or real values; when the "pcf" is a density estimate
   based on data "dendat", then specifying a value for the parameter "dendat",
   we can plot the data together with the estimate of a level set}
  \item{xaxt}{a character which specifies the x axis type; either "s" or "n";
   see "par"}
  \item{yaxt}{a character which specifies the y axis type; either "s" or "n";
   see "par"}
  \item{cex.axis}{the magnification to be used for axis annotation}
}
%\details{}
\value{
Makes a plot to the graphics window
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemela }
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
