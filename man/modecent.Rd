\name{modecent}
\alias{modecent}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ 
Returns locations of modes of a density estimate  
}
\description{
With "modecent" we may find the location of the modes of a function.
Argument to function "modecent" is the level set tree of the function.
A mode is defined here to be the barycenter of a set at which
the estimate has a local maximum. 
}
\usage{
modecent(lst)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{lst}{level set tree object;
  list containing lst$parent and lst$center;
  output of "leafsfirst", "profgene", "profhist", or "profkern"  }
}
%\details{  }
\value{ modenum*d-matrix; for each mode its location.
  The rows of output are d-vectors specifying the locations of the modes.}
%\references{ http://www.rni.helsinki.fi/~jsk/denpro}
\author{ Jussi Klemelä }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ 
\code{\link{leafsfirst}},
\code{\link{profkern}}, 
\code{\link{plotbary}} 
}

\examples{
set.seed(1)
dendat<-matrix(rnorm(20),10)   #10*2 data-matrix
pcf<-pcf.kern(dendat,h=2,N=c(10,10))
lst<-leafsfirst(pcf)
modecent(lst)

#set.seed(1)
#dendat<-matrix(rnorm(20),10)   #10*2 data-matrix
#ph<-profhist(dendat,binlkm=3)
#modecent(ph)           

#          [,1]      [,2]
#[1,] 2.0004324 -1.593620
#[2,] 0.3798261  0.890701

#recnum<-3
#d<-2
#values<-seq(1:recnum)
#recs<-matrix(0,recnum,2*d)
#recs[1,]<-c(0,1,0,1)
#recs[2,]<-c(0,1,1,2)
#recs[3,]<-c(1,2,0,2)
#pg<-profgene(values,recs)    
#modecent(pg)

#     [,1] [,2]
#[1,]  1.5    1  

}
\keyword{multivariate}%-- one or more ...
