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
profhist(dendat, binlkm, cvol=T, ccen=T, cfre=F)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{dendat}{n*d data matrix}
 \item{binlkm}{positive integer;
  number of bins of the histogram estimate in one direction. 
  We use the same number of bins for every direction, thus the total
  number of bins is binlkm to the power of d}
 \item{cvol}{T or F; TRUE if one wants that the volumes of
 the separated parts of the level sets are returned, 
 note that one needs volumes for drawing the volume plot}
 \item{ccen}{T or F; TRUE if one wants that the barycenters of the 
 separated parts of the level sets are returned}
 \item{cfre}{T or F; TRUE if one wants that the frequencies of the
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
\author{ Jussi Klemelä }

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
#ph

#$parent
#     [,1] [,2] [,3] [,4]
#[1,]    0    0    1    3
#
#$level
#           [,1]       [,2]      [,3]      [,4]
#[1,] 0.09935159 0.09935159 0.1987032 0.3974064
#
#$invalue
#           [,1]       [,2]      [,3]      [,4]
#[1,] 0.09935159 0.09935159 0.1987032 0.3974064
#
#$volume
#         [,1]     [,2]     [,3]     [,4]
#[1,] 5.032632 1.006526 2.013053 1.006526
#
#$center
#           [,1]      [,2]        [,3]      [,4]
#[1,] -0.1063558  2.000432 -0.02532547 0.3798261
#[2,]  0.6422689 -1.593620  0.26962082 0.8907010
#
#$nodefrek
#          [,1]       [,2]      [,3]      [,4]
#[1,] 0.8941643 0.09935159 0.5961095 0.3974064
#
#$recs
#            [,1]        [,2]       [,3]       [,4]
#[1,] -0.83562861 -0.02532547  1.5117812  2.7539415
#[2,] -0.02532547  0.78497766  0.2696208  1.5117812
#[3,] -0.83562861 -0.02532547 -0.9725395  0.2696208
#[4,]  1.59528080  2.40558394 -2.2146999 -0.9725395
#[5,] -0.02532547  0.78497766 -0.9725395  0.2696208
#[6,] -0.83562861 -0.02532547  0.2696208  1.5117812
#
#$hisfrekv
#           [,1]      [,2]      [,3]       [,4]       [,5]       [,6]
#[1,] 0.09935159 0.3974064 0.1987032 0.09935159 0.09935159 0.09935159
#
#$lsets
#     [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]    1    0    1    1    1    1
#[2,]    0    1    0    0    0    0
#[3,]    0    0    0    0    1    1
#[4,]    0    0    0    0    0    1

}

\keyword{smooth}% at least one, from doc/KEYWORDS
%\keyword{ Graphics}% __ONLY ONE__ keyword per line





