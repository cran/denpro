\name{plotmodet}
\alias{plotmodet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Plots a mode graph  
}
\description{
Plots one window of a mode graph.
}
\usage{
plotmodet(mt, coordi = 1, colot = NULL, shift = 0, xlim = NULL, 
xlab = "", ylab = "", horilines = NULL, symbo = 20, loga = NULL, 
lty = "dashed", cex.axis = 1, title = TRUE, cex.sub = 1, cex.lab = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mt}{mode graph; output of "modegraph"}
  \item{coordi}{integer 1,...,d }
  \item{colot}{colors }
  \item{shift}{real number; to shift the parent-child lines in the case of 
   overlays }
  \item{xlim}{2-vector of reals, giving the interval of x-axis  }
  \item{xlab}{ labeling for the x-axis}
  \item{ylab}{ labeling for y-axis }
  \item{horilines}{ a vector of indexes; horizonta lines will be plotted }
  \item{symbo}{ labeling for modes; "L" or "N" }
  \item{loga}{ "y" or "n"; whether logarithmic y-axis will be used }
  \item{lty}{ line type }
  \item{cex.axis}{magnification factor for the axis annotation; see "par" }
  \item{title}{TRUE or FALSE; makes a subtitle which shows the 
   value of "coordi" }
  \item{cex.sub}{magnification factor for the subtitle; see "par" }
  \item{cex.lab}{magnification factor for the axis labels (names of x and y); 
   see "par" }
}
%\details{}
\value{
Plots to the graphics device
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }

\seealso{
\code{\link{modegraph}}, 
\code{\link{lstseq.kern}}, 
}
\examples{
dendat<-sim.data(n=200,type="mulmod")

h1<-0.9
h2<-2.2
lkm<-5
hseq<-hgrid(h1,h2,lkm)

N<-c(16,16)
estiseq<-lstseq.kern(dendat,hseq,N,lstree=TRUE)

mt<-modegraph(estiseq)

plotmodet(mt)

plotmodet(mt,coordi=2)

}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{hplot }% __ONLY ONE__ keyword per line
