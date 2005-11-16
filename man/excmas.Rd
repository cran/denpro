\name{excmas}
\alias{excmas}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Calculates the excess masses associated with the 
nodes of a level set tree
}
\description{
Given a level set tree, 
calculates the excess mass for every node of the tree. 
The excess mass is the probability mass under the curve, 
over a given level, 
in a certain branch of the level set tree.
For example, function "leafsfirst" returns a level set tree. 
}
\usage{
excmas(lst)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{lst}{level set tree object;
           list of vectors which is the output of for example 
           "leafsfirst", "profkern", "profhist", or "profgene" }
}
%\details{}
\value{
  Vector whose length is equal to the number of nodes of the level set tree.
}
%\references{{http://www.rni.helsinki.fi/~jsk/denpro}
\author{ Jussi Klemelä }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ 
\code{\link{leafsfirst}},
\code{\link{profkern}}, 
\code{\link{plottree}}
}

\examples{
set.seed(1)
dendat<-matrix(rnorm(20),10)   #10*2 data-matrix
pcf<-pcf.kern(dendat,h=2,N=c(10,10))
lst<-leafsfirst(pcf)
excmas(lst)

#set.seed(1)
#dendat<-matrix(rnorm(20),10)   #10*2 data-matrix
#pk<-profkern(dendat,h=1,N=c(8,8),Q=4)
#excmas(pk)

#     [,1]      [,2]      [,3]       [,4]
#[1,]    1 0.9833526 0.3511974 0.07023947

#set.seed(1)
#dendat<-matrix(rnorm(20),10)   #10*2 data-matrix
#ph<-profhist(dendat,binlkm=3)
#excmas(ph)

#     [,1] [,2] [,3] [,4]
#[1,]  0.9  0.1  0.4  0.2

}
\keyword{multivariate}%-- one or more ...
