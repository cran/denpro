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
pcf<-pcf.kern(dendat,h=2,N=c(16,16))
lst<-leafsfirst(pcf)
modecent(lst)

}
\keyword{multivariate}%-- one or more ...
