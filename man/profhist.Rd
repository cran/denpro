\name{profhist}
\alias{profhist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Calculates the level set tree of a histogram
}
\description{
Given a data matrix, 
returns the level set tree of a histogram, 
which is constructed from bins of equal size.
}
\usage{
profhist(dendat, binlkm, cvol=TRUE, ccen=TRUE, cfre=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{dendat}{n*d data matrix}
 \item{binlkm}{positive integer;
  number of bins of the histogram estimate in one direction. 
  We use the same number of bins for every direction, thus the total
  number of bins is binlkm to the power of d}
 \item{cvol}{TRUE if one wants that the volumes of
 the separated parts of the level sets are returned, 
 note that one needs volumes for drawing the volume plot}
 \item{ccen}{TRUE if one wants that the barycenters of the 
 separated parts of the level sets are returned}
 \item{cfre}{TRUE if one wants that the frequencies of the
 separated parts of the level sets are returned }
}
%\details{ }
\value{
  An augmented level set tree.  The level set tree is a list of
  vectors.  The elements of the vectors supply information for each
  node of the tree.  Below we denote with "nodenum" the number of
  nodes of the tree.
  \item{parent}{"nodenum"-vector of integers in range 0,..., nodenum-1;
  links to the parent of each node. Root nodes are marked with 0.}
  \item{level}{"nodenum"-vector of positive real numbers;
  level of the level set from which the set corresponding to the node
  is a part of.}
  \item{invalue}{"nodenum"-vector of positive integers;
    level of the level set in terms of the original frequencies 
    (these values are not normalized so that estimate would integrate to one}
  \item{volume}{"nodenum"-vector of positive real numbers; 
   gives volumes for sets corresponding to each node}
  \item{center}{d*nodenum-matrix; gives the barycenter for sets 
    corresponding to each node}
  \item{nodefrek}{"nodenum"-vector of positive integers;
    number of observations in the set corresponding to the node. 
    This is useful in cluster analysis applications.}
  \item{recs}{atomnum*(2*d)-matrix, where atomnum is the number of
   non-empty bins. Matrix defines the non-empty bins,
   these are the "atoms" of level sets}
  \item{hisfrek}{atomnum-vector of  positive integers;
     number of observations in non-empty bins}
  \item{lsets}{nodenum*atomnum-matrix: describes the sets associated
  with the nodes. We have 1 in column "c" if the atom described in c:th 
  row of recs is a part of the set associated with this node }
}
%\references{http://www.rni.helsinki.fi/~jsk/denpro}
\author{ Jussi Klemela }

\note{ Applies the naive algorithm of pairwise comparison of the 
separated components of the level sets, to find which components 
touch each other.}

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ 
\code{\link{plotvolu}},
\code{\link{plotbary}}, 
\code{\link{plottree}}, 
\code{\link{profgene}} 
}

\examples{
set.seed(1)
dendat<-matrix(rnorm(20),10)   #10*2 data-matrix
ph<-profhist(dendat,binlkm=3,cfre=TRUE)

}

\keyword{smooth}% at least one, from doc/KEYWORDS
%\keyword{ Graphics}% __ONLY ONE__ keyword per line





