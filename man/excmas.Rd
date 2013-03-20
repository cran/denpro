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
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ 
\code{\link{leafsfirst}},
\code{\link{plottree}}
}

\examples{
set.seed(1)
dendat<-matrix(rnorm(20),10)   #10*2 data-matrix
pcf<-pcf.kern(dendat,h=1,N=c(16,16))
lst<-leafsfirst(pcf)
excmas(lst)

}
\keyword{multivariate}%-- one or more ...
\keyword{smooth}
