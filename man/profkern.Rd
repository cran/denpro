\name{profkern}
\alias{profkern}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Calculates a level set tree of a kernel estimate
}
\description{
Given a data matrix, returns a level set tree of a kernel estimate, 
when the dimension of the estimate is less or equal to 4.
}
\usage{
profkern(dendat, h, N, Q, cvol=TRUE, ccen=TRUE, cfre=FALSE,
kernel="epane", compoinfo=FALSE, trunc=3, threshold=0.0000001,
sorsa="crc",hw=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{dendat}{n*d data matrix; d<=4}
 \item{h}{positive real number or vector of positive real numbers; 
      smoothing parameters from the largest to the smallest}
 \item{N}{d-vector of dyadic integers >=4;
          number of knots for every direction: kernel
          estimate will be evaluated on these knots}
 \item{Q}{positive integer; number of quantization levels} 
 \item{cvol}{TRUE or FALSE; TRUE if one wants that the volumes of
 the separated parts of the level sets are returned, 
 note that one needs the volumes for drawing a volume plot}
 \item{ccen}{T or F; TRUE if one wants that the barycenters of the 
 separated parts of the level sets are returned,
 note that one needs the barycenters for drawing a barycenter plot}
 \item{cfre}{T or F; TRUE if one wants that the frequencies for 
 separated parts of the level sets are returned (not implemented)}
 \item{kernel}{either "epane" or "gauss"; 
  the kernel of the kernel estimate}
 \item{compoinfo}{TRUE or FALSE; whether definitions of separated
  regions of level sets will be returned}
 \item{trunc}{positive real number; truncation of the Gaussian kernel
  to the interval [-trunc,trunc]}
 \item{threshold}{positive real number; kernel will be truncated to
  0 if the value is smaller than the threshold (not implemented)}
 \item{sorsa}{if not "crc" uses slower R code}
 \item{hw}{lower level parameter}
}
%\details{}
\value{
  Returns a list of level set trees: for each value  of h
  the level set tree of the corresponding kernel estimate. 
  If h is not a vector but scalar, then only one
  level set tree will be returned.
  The level set tree is a list of vectors. 
  The elements of the vectors supply information for each node of the tree. 
  Below we denote with "nodenum" the number of nodes of the tree.
  \item{parent}{"nodenum"-vector of integers in the range 0,..., nodenum-1;
  links to the parent of each node. Root nodes are marked with 0.}
  \item{level}{"nodenum"-vector of positive real numbers;
  level of the level set from which the set corresponding to the node
  is a part of.}
  \item{volume}{"nodenum"-vector of positive real numbers; 
   volumes of sets corresponding to each node}
  \item{center}{d*nodenum-matrix; barycenters of sets 
    corresponding to each node}
  \item{invalue}{"nodenum"-vector of positive integers;
    level of the level set in terms of original frequencies 
    (these values are not normalized so that the 
    estimate would integrate to one}
  \item{component}{Will be returned if "compoinfo"=TRUE; for each
  node of the level set a tree pointer to "AtomlistAtom"}
  \item{AtomlistAtom}{Will be returned if "compoinfo"=TRUE; pointers
  to "index"}
  \item{AtomlistNext}{Will be returned if "compoinfo"=TRUE; pointers
  which define the list of atoms for each separated component (for each node
  of the level set tree)}
  \item{index}{Will be returned if "compoinfo"=TRUE; matrix with d
  columns: determines a knot of the multivariate grid where estimate
  was evaluated}
  \item{nodefrek}{"nodenum"-vector of positive integers;
    number of observations in the set corresponding to the node. 
    This is useful in cluster analysis applications. (Not implemented.)}
}
%\references{http://www.rni.helsinki.fi/~jsk/denpro}
\author{ Jussi Klemelä }

\note{Applies the DynaDecompose algorithm, described in
the article
"Algorithms for manipulation of level sets of nonparametric density
estimates",
by Jussi Klemelä }

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ 
\code{\link{plotvolu}},
\code{\link{plotbary}}, 
\code{\link{plottree}}, 
}

\examples{
set.seed(1)
dendat<-matrix(rnorm(20),10)      # 10*2 data-matrix
pk<-profkern(dendat,h=1,N=c(8,8),Q=4)
plotvolu(pk)

dendat<-sim.data(n=200,type="mulmod")
pk<-profkern(dendat,h=1,N=c(64,64),Q=30)
plotvolu(pk)

}

\keyword{smooth}% at least one, from doc/KEYWORDS
















