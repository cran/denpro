\name{scaletable}
\alias{scaletable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implements the scale and shape table
}
\description{
Implements the scale and shape table, which is a dynamic tool to visualize a
scale of multivariate estimates and to inspect the shape of the estimates.
Takes as an input a scale of estimates. 
}
\usage{
scaletable(estiseq, paletti = NULL, shift = 0, ptext = 0, 
ptextst = 0, bm = NULL, levnum = 60, levnumst = 60, 
redu=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estiseq}{a list of density estimates, contains 5 fields: 
                 pcfseq, lstseq, stseq, hseq, and type}
  \item{paletti}{a character vector of color names}
  \item{shift}{a real number; possible shifting of mode tree lines}
  \item{ptext}{a real number; lifting of the mode labels}
  \item{ptextst}{a real number; lifting of the tail labels}
  \item{bm}{branching map; just to fasten the plotting in the case the
  branching map has been already calculated}
  \item{levnum}{positive ineteger}{the number of levels in level set trees}
  \item{levnumst}{positive ineteger}{the number of levels in shape trees}
  \item{redu}{T or F; whether the number of levels of levels will be reduced,
   see argument "levnum"}
}
%\details{}
\value{
Plots 7 frames: 
control frame,
map of branches and mode graph (the biggest frames), 
volume plot, and barycenter plot (2 medium frames).
radius plot, and location plot (2 small frames).

One chooses from the control panel the frame which one wants to manipulate.
To return to the control panel, 
click at the bottom of the frame which is currently active.

Interaction in the mode graph frame:
At the beginning the mode tree shows the first coordinate.
Click at the top of the mode tree frame to change the coordinate
of the mode tree.
Click with the mouse at the scale of smoothing parameter values;
the frames of 
volume plot barycenter plot, radius plot, and location plot
change to correspond to this smoothing parameter value.
The vertical position of the click point determines the smoothing parameter
value, and the horizontal position does not matter.
(At the beginning the frames correspond
to the estimate with the largest smoothing parameter,
and radius plot and location plot show the 10\% level set.)

Interaction in the volume plot frame:
Click inside of the figure to zoom in (to change the range of x-coordinates).
Click at the top of the frame to return to the original scale.
Click at the left to the y-coordinate axis to choose the level of the 
radius plot and location plot.

Interaction in the frame of the barycenter plot:
Click at the top of the frame or inside the frame to change the coordinate.
(At the beginning the barycenter plot shows the first coordinate.)

Interaction in the frame of the radius plot:
Choose with the mouse click at the top of the frame the reference point
of the radius and location plot.
The reference point is either th elocation of the maximum of the estimate,
or the barycenter of the level set.

Interaction in the frame of the location plot:
Click at the top of the frame or inside the frame to change the coordinate.
(At the beginning the location plot shows the first coordinate.)

Interaction in the frame of the map of branches:
Rotate the perspective plot with the mouse clicks.
Each click rotates 10 degrees.
The direction of the rotation is up or down, left or right, dpending on the 
position of the mouse click relative to the center of the frame.

To finish the program, click with the mouse at the label "STOP"
at the control frame.
}

\author{Jussi Klemel\"a}
%\note{}

\seealso{\code{\link{lstseq.kern}}}

\examples{
dendat<-sim.data(n=200,type="mulmod")

h1<-0.9
h2<-2.2
lkm<-5
hseq<-hgrid(h1,h2,lkm)

N<-c(16,16)
estiseq<-lstseq.kern(dendat,hseq,N,lstree=TRUE,level=0.1) 

# scaletable(estiseq)

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
