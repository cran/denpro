\name{locofmax}
\alias{locofmax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the location of the maximum of a function
}
\description{
Calculates the location of the maximum of a piecewise constant function.
}
\usage{
locofmax(pcf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{piecewise constant function; output of "pcf.kern" or "eval.func"}
}
%\details{}
\value{
a d-vector; the location of the maximum;
this is defined as the barycenter of the set where the function has the maximum
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{ 
\code{\link{pcf.kern}},
\code{\link{eval.func}}
}
\examples{
dendat<-sim.exa(n=100)
pcf<-pcf.kern(dendat,h=1,N=c(25,25))

locofmax(pcf)

}
\keyword{multivariate }% at least one, from doc/KEYWORDS

