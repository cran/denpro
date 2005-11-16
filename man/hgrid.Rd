\name{hgrid}
\alias{hgrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns a grid of smoothing parameter values }
\description{
Returns a grid of smoothing parameter values for kernel estimates,
either with a logarithmic spacing, or with equal spacing,
in decreasing order.
}
\usage{
hgrid(h1, h2, lkm, base = 10)
}

\arguments{
  \item{h1}{ the lowest smoothing parameter value }
  \item{h2}{ the largest smoothing parameter value }
  \item{lkm}{ the number of smoothing parameters in the grid }
  \item{base}{ the base of the logarithm, used in the logarithmic spacing  }
}

\value{
a vector of smoothing parameters
}
\author{ Jussi Klemelä }

\seealso{ 
\code{\link{lstseq.kern}},
\code{\link{branchmap}},
\code{\link{modegraph}}, 
\code{\link{exmap}}
}

\examples{
h1<-0.9
h2<-2.2
lkm<-5
hseq<-hgrid(h1,h2,lkm)

}
\keyword{ smooth } % at least one, from doc/KEYWORDS


