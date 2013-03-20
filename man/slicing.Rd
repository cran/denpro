\name{slicing}
\alias{slicing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns a one- or two-dimensional slice of a multivariate function }
\description{
Returns a one- or two-dimensional slice of a multivariate function.
The functions are represented as piecewise constant functions.
}
\usage{
slicing(pcf, vecci, d1 = 1, d2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcf}{ a piecewise cosntant function }
  \item{vecci}{ vector inside the support of the function
   which fixes the variables;
   length of "vecci" is d-2 or d-1, depending on whether the slice is
   2- or 1-dimensional}
  \item{d1}{ integer 1,...,d; the first variable of the slice }
  \item{d2}{ integer 1,...,d, not equal to d1; the second variable of the 
   slice }
}
%\details{}
\value{
a piecewise constant function
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemela}
%\note{ ~~further notes~~ }
\seealso{ 
\code{\link{draw.pcf}} 
}
\examples{
# 2D slice of a 3D function

N<-c(26,26,26)  # choose the grid size
copula<-"gauss"
margin<-"student"
b<-4
support<-c(-b,b,-b,b,-b,b)
r<-0.5       # parameter of the copula
t<-c(2,2,3)  # degreeds of freedom for the student margin
pcf<-pcf.func(copula,N,marginal=margin,support=support,r=r,t=t)

sl<-slicing(pcf,d1=1,d2=2,vecci=0)
dp<-draw.pcf(sl)
persp(dp$x,dp$y,dp$z,theta=30,phi=30)   

# 1D slice of a 2D function
N<-c(60,60)                  
pcf<-sim.data(N=N,type="mulmod")

sl<-slicing(pcf,vecci=0)
dp<-draw.pcf(sl)
plot(dp$x,dp$y,type="l")

}

\keyword{ multivariate }



