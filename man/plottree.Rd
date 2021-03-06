\name{plottree}
\alias{plottree}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Makes a tree plot of a level set tree,
of a shape tree, or of a tail tree
}
\description{
Plots a tree plot of a level set tree, of a shape tree, 
or of a tail tree to the graphics window.
}
\usage{
plottree(lst, 
plot=T, data=F, crit=NULL, orderrule="distcenter",
modelabel=TRUE, ptext=0, leimat=NULL, symbo=NULL, 
info=NULL, infolift=0, infopos=0, infochar=NULL,
xmarginleft=0, xmarginright=0, ymargin=0, 
xlim=NULL, ylim=NULL, col="black", col.axis="black", 
linecol=rep("black",length(lst$parent)), 
pch=21, dimen=NULL, yaxt="s", axes=T,
cex=NULL, nodemag=NULL, linemag=1, 
cex.axis = 1, ylab = "", cex.lab = 1, colo = FALSE, paletti = NULL,
lowest = "dens")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{lst}{level set tree; list of vectors.
   The list contains at least vectors "level", "volume", and "parent". 
   For example, functions "leafsfirst", "profkern" and "profhist" 
   give a level set tree as an output.
  }
  \item{plot}{T or F; TRUE if we make a plot, otherwise FALSE.}
  \item{data}{T or F; TRUE if we want output to contain some
   information, for example an ordering for siblings. 
   This option is needed only by other plotting
   functions of the package, it is not needed by the end user.}
  \item{crit}{d-vector of real numbers; gives a way to control
   ordering of siblings. The leftmost sibling is the one whose
   barycenter is furthest away from vector "crit", 
   in the Euclidean distance.}
  \item{orderrule}{lower level paremeter}
  \item{modelabel}{T or F; TRUE if the modes will be labelled.
   The default is to use labels M1, M2,...}
  \item{ptext}{non-negative real number; amount by which the mode
   labels will be lifted.}
  \item{leimat}{vector of characters; length of the vector should be
   equal to the number of modes of the estimate.
   This option is for the case we do not want the ordering of the
   labels to be done automatically.}
  \item{symbo}{character; for example "L". The default value for the
   automatic labelling of the modes is to use M1, M2,...
   With "symbo" we may switch to L1, L2,..., for example.}
  \item{info}{vector of numbers or characters, whose length is equal
   to the number of nodes of the level set tree.
   The elements of "info" will be placed on the right side of the 
   nodes.
   For example "info" may be generated by "excmas" or we may define "info"
   to contain the frequencies of the nodes.
   (Frequencies may be obtained directly from the function "profhist".)}
  \item{infolift}{real number; controls the vertical positioning of
   the elements of "info".}
  \item{infopos}{real number; controls the horizontal positioning of
   the elements of "info". Negative "infopos" will move elements of
   "info" to the right hand side.}
   \item{infochar}{vector of characters whose length is equal to the
   number of nodes of the level set tree; we may annotate nodes with
   some infromation}
  \item{xmarginleft}{nonnegative real number; adds more margin on the 
   left hand side. The box around the plot will be moved to the left
   with the amount "xmarginleft".}
  \item{xmarginright}{nonnegative real number; adds more margin on the 
   right hand side. The box around the plot will be moved to the right
   with the amount "xmarginright".}
  \item{ymargin}{nonnegative real number; adds more margin on the 
   top of the plot. The box around the plot will be moved up
   with the amount "ymargin".}
  \item{xlim}{vector of 2 real numbers; gives the limits for the scale
   of x-axis.}
  \item{ylim}{vector of 2 real numbers; gives the limits for the scale
   of y-axis.}
  \item{col}{colour for the nodes; for example "black" or "blue".}
  \item{col.axis}{colour for the x and y-axis; for example "black" or "blue".}
  \item{linecol}{colour for the lines joining nodes; 
  for example "black" or "blue".}
  \item{pch}{symbol for the nodes of the tree; integer 19-25; 
   see help(points) for the definitions}
  \item{dimen}{positive integer; number of dimensions of the estimate.
   If "dimen" is not NULL, then we plot without using the usual
   ordering of the siblings. This option is used when we do not have
   barycenters available.}
  \item{yaxt}{axis type, see "par"}
  \item{axes}{T or F}
  \item{cex}{ magnification factor for the symbols; see "par" }
  \item{nodemag}{ magnification factor for the points }
  \item{linemag}{ magnification factor for the lines joining the points }
  \item{cex.axis}{magnification factor for the axis annotation; see "par" }
  \item{ylab}{character; y-axis label (name of the y-variable); 
   for example in a tree plot of a level set tree this could be "level"}
  \item{cex.lab}{magnification factor for the axis labels (names of x and y); 
   see "par" }
  \item{colo}{TRUE or FALSE; if TRUE, nodes and lines joining the nodes will
   be colored}
  \item{paletti}{a vector of names of colors; these colors will be used to
   color the nodes and the lines joining the nodes}
  \item{lowest}{character string; if lowest="dens", then y-axis starts at 0}  
}
\value{
  By default a plot is made on the graphics window. 
}

%\references{http://www.rni.helsinki.fi/~jsk/denpro}
\author{ Jussi Klemela }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ 
\code{\link{plotvolu}}, 
\code{\link{plotbary}},
\code{\link{leafsfirst}}  
}

\examples{
dendat<-sim.data(n=50,type="mulmod")

pcf<-pcf.kern(dendat,h=1,N=c(32,32))
lst<-leafsfirst(pcf)
td<-treedisc(lst,pcf,ngrid=30)
plottree(td)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line





