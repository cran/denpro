\name{colors2data}
\alias{colors2data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assigns colors to data points
}
\description{
Assigns colors to data points according to cluster membership, 
when the clusters are defined as high density regions 
}
\usage{
colors2data(dendat, pcf, lst, paletti = NULL, clusterlevel = NULL, nodes = NULL,
type = "regular")
}
\arguments{
  \item{dendat}{n*d matrix of real numbers; the data matrix.}
  \item{pcf}{piecewise constant function; output of "pcf.kern", or
  "pcf.func", for example.}
  \item{lst}{level set tree; a list of vectors.
   The list contains at least vectors "level", "volume", and "parent". 
   For example, function "leafsfirst" computes a level set tree.}
  \item{paletti}{a vector of names of colors; for example paletti=c("black","red","blue").}
  \item{clusterlevel}{a positive real number; gives the level which is such that the
   disconnected components of the level set with this level define the clusters.}
  \item{nodes}{a vector of positive integers; contains pointers to the nodes of
   level set tree "lst"; the nodes indicate which disconnected components of level sets
   define the clusters.}
  \item{type}{either "regular" or "adaptive"; if type="regular", then the piecewise
   constant function has a regular partition; if type="adaptive", then the piecewise
   constant function has an adaptive partition.} 
}
%\details{}
\value{
a list with components "datacolo" and "ord".
"datacolo" is a vector of color names.
This vector assigns a color to each row of matrix "dendat".
"ord" is a vector of positive integers.
This vector contains pointers to rows of the matrix "dendat".
The pointers indicate the order so that the data points with the lowest estimated
density are first.
}
%\references{http://jussi.klemela.com/denpro/}
\author{Jussi Klemela}
%\note{ ~~further notes~~ }

 \seealso{
\code{\link{pcf.kern}}
\code{\link{leafsfirst}}
}

\examples{
# generate data
seed<-1
n<-50
d<-2
l<-3; D<-4; c<-D/sqrt(2)
M<-matrix(0,l,d); M[2,]<-c; M[3,]<--c
sig<-matrix(1,l,d)
p<-rep(1/l,l)
dendat<-sim.data(type="mixt",n=n,M=M,sig=sig,p=p,seed=seed)

# colored volume function 
h<-(4/(d+2))^(1/(d+4))*n^(-1/(d+4))*apply(dendat,2,sd)
N<-rep(20,d)
pcf<-pcf.kern(dendat,h,N)
lst<-leafsfirst(pcf)

# colored scatter plot
cd<-colors2data(dendat,pcf,lst)
plot(dendat,col=cd$datacolo)

# colored scatter plot: plot the observations so that the observations with
# the highest density value are plotted last
cd<-colors2data(dendat,pcf,lst)
plot(dendat[cd$ord,1],dendat[cd$ord,2],col=cd$datacolo[cd$ord])

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line


