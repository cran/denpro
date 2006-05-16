\name{lstseq.kern}
\alias{lstseq.kern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates a scale of kernel estimates
}
\description{
Calculates a scale of kernel estimates corresponding to a scale of smoothing 
parmeters. 
}
\usage{
lstseq.kern(dendat, hseq, N, lstree = NULL, level = NULL,
Q = NULL, kernel = "epane", hw = NULL, algo = "leafsfirst")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendat}{n*d matrix of real numbers; the data matrix}
  \item{hseq}{a vector of positive real numbers;
   the sequence should be monotonic}
  \item{N}{vector of d positive integers; the dimension of the grid where the
    kernel estimate will be evaluated; we evaluate the estimate on a regular
    grid which contains the support of the kernel estimate}
  \item{lstree }{ if NULL, then level set trees are not calculated }
  \item{level}{ NULL or a real number between 0 and 1; 
   if NULL, then shape trees are not calculated;
   if number, then it is the level in percents of the maximum
   of the level sets for which the shape trees are calculated }
  \item{Q}{positive integer; needed only in the DynaDecompose algorithm, see
   parameter "algo"; the number of levels in the level set trees}
  \item{kernel}{"epane" or "gauss"; the kernel is either the
   Bartlett-Epanechnikov product kernel or the standard Gaussian}
  \item{hw}{positive integer; parameter for time localized kernel estimation;
   gives the smoothing parameter for the temporal smoothing}
  \item{algo}{"leafsfirst" or "dynadecompose"}
}
%\details{}
\value{
  A list with components
  \item{lstseq}{a list of level set trees}
  \item{pcfseq}{a list of piecewise constant functions}
  \item{stseq}{a list of shape trees}
  \item{hseq}{a vector of smoothing parameters corresponding to the members
  in the sequences}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jussi Klemelä}
%\note{ ~~further notes~~ }

\seealso{
\code{\link{scaletable}}
}
\examples{
dendat<-sim.data(n=200,type="mulmod")

h1<-0.9
h2<-2.2
lkm<-5
hseq<-hgrid(h1,h2,lkm)

N<-c(16,16)
estiseq<-lstseq.kern(dendat,hseq,N)

}
\keyword{smooth}% at least one, from doc/KEYWORDS

