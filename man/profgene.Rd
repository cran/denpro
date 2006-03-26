\name{profgene}
\alias{profgene}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ 
Calculates the level set tree of a rectangularwise constant function
}
\description{
Returns the level set tree of a function
which is rectangularwise constant.
}
\usage{
profgene(values, recs, frekv=NULL, cvol=TRUE, ccen=TRUE, cfre=FALSE,
outlsets=TRUE, invalue=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{values}{recnum-vector of positive real values; 
 we want to define a rectangularwise
 constant non-negative function and we denote with "recnum" the number of
 rectangles at which the function has values greater than 0.
 }
 \item{recs}{recnum*(2*d)-matrix; corresponding to each element of vector
 "value" we specify a rectangle. 
 In each row of the "recs" we specify
 a rectangle at which the estimate has the value given at the corresponding
 element of "value". 
 Rows of "recs" have the form c(b1,e1,...,bd,ed)
 when rectangle is the product of intervals [bi,ei], i=1,...,d. }
 \item{frekv}{recnum-vector; for each rectangle the number of observations
 in this rectangle. 
 Supplying this argument is useful in clustering: 
 one may see the number of observations in each node of the level set tree. }
 \item{cvol}{TRUE if one wants that the volumes of
 separated parts of the level sets are returned,
 note that one needs volumes for drawing volume plots }
 \item{ccen}{TRUE if one wants that the barycenters of the
 separated parts of the level sets are returned }
 \item{cfre}{TRUE if one wants the frequencies of
 separated parts of the level sets to be returned }
 \item{outlsets}{not needed}
 \item{invalue}{not needed}
}
%\details{}
\value{
  An augmented level set tree. 
  The level set tree is a list of vectors. 
  The elements of the vectors supply information for each node of the tree. 
  Below we denote with "nodenum" the number of nodes of the tree.
  \item{parent}{"nodenum"-vector of integers in range 0,..., nodenum-1;
  links to the parent of each node. Root nodes are marked with 0.}
  \item{level}{"nodenum"-vector of positive real numbers;
  level of the level set from which the set corresponding to the node
  is a part of.}
  \item{volume}{"nodenum"-vector of positive real numbers; 
   volume of sets corresponding to each node}
  \item{center}{d*nodenum-matrix; barycenters of sets 
    corresponding to each node}
  \item{invalue}{"nodenum"-vector of positive integers;
    level of the level set in terms of original frequencies 
    (these values are not normalized so that estimate would integrate to one)}
  \item{nodefrek}{"nodenum"-vector of positive integers;
    number of observations in the set corresponding to node. 
    Useful in cluster analysis applications.}
  \item{lsets}{nodenum*binnum-matrix; describes the sets associated
  with nodes. We have 1 in column "c" if the bin described in c:th 
  row of recs is part of the set associated with this node }
}
%\references{http://www.rni.helsinki.fi/~jsk/denpro}
\author{ Jussi Klemelä }
\note{ Applies the naive algorithm of pairwise comparison of the 
separated components of the level sets, to find which components 
touch each other. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{  
\code{\link{profhist}},\code{\link{profkern}},\code{\link{plotvolu}} 
}

\examples{
recnum<-3
d<-2
value<-seq(1:recnum)
recs<-matrix(0,recnum,2*d)
recs[1,]<-c(0,1,0,1)
recs[2,]<-c(0,1,1,2)
recs[3,]<-c(1,2,0,2)
pg<-profgene(value,recs)

}
\keyword{multivariate}%-- one or more ...
