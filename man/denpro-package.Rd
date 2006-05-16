\name{denpro-package}
\alias{denpro-package}
\alias{denpro}
\docType{package}
\title{
Visualization of multivariate, functions, sets, and data
}
\description{
The package provides tools to (1) visualize multivariate
density functions and density estimates with level set trees,
(2) visualize level sets with shape trees, (3) visualize
multivariate data with tail trees, (4) visualize scales of
multivariate density estimates with mode graphs and branching
maps, and (5) visualize anisotropic spread with 2D volume
functions and 2D probability conetent functions. With level set
trees one visualizes mode structure, with shape trees one
visualizes shapes of level sets of unimodal densities, and with
tail trees one visualizes connected data sets. The kernel
estimator is implemented but the package may be applied for
visualizing other density estimates, which have to be setwise
constant.
}
\details{
\tabular{ll}{
Package: \tab denpro\cr
Version: \tab 0.8.2\cr
Date: \tab 2006-05-16\cr
Depends: \tab R\cr
License: \tab GPL version 2 or newer\cr
URL: \tab http://www.rni.helsinki.fi/~jsk/denpro\cr
Packaged: \tab Thu May 18 15:34:27 2006; jsk\cr
Built: \tab R 2.2.1; i686-pc-linux-gnu; 2006-05-18 15:34:43; unix\cr
}


The main function:

    *leafsfirst 

Functions to plot a level set tree, shape tree, or tail tree:

    *plotvolu 
    *plotbary
    *plottree 
    *treedisc 
    *prunemodes 

Visualization of scales of estimates:

    * lstseq.kern 
    * modegraph 
    * plotmodet 
    * branchmap 
    * plotbranchmap 
    * scaletable 
    * exmap 
    * plotexmap 
    * hgrid :

Visualization of anisotropic spread 
(2D volume function and 2D probability content function):

    * stseq 
    * shape2d 
    * plotvolu2d 
    * plotdelineator

Other utilities:

    * pcf.kern 
    * tree.segme
    * lst2xy
    * graph.matrix
    * plotbary.slide
    * Calculation of level set trees with "DynaDecompose" algorithm 
      and with naive algorithm
          o profkern 
          o profhist 
          o profgene 
          o profeval
    * locofmax 
    * modecent 
    * excmas 

Other visualization tools:

    * draw.pcf 
    * slicing 
    * draw.levset 
    * paracoor 
    * pp.plot
    * qq.plot
    * plot.kernscale
    * dist.func
    * paraclus

Miscallenous:

    * pcf.func 
    * sim.data 
    * dend2parent
    * explo.compa
    * nn.radit
    * nn.likeset
    * rotation
    * liketree
    * leafsfirst.visu
    * plottwin

Index (alphaphetical):
\preformatted{
branchmap               Calculates a branching map from a sequence of
                        level set trees
draw.levset             Plots a level set of a 2D function
draw.pcf                Prepares the plotting of a 2D or 1D piecewise
                        constant function
excmas                  Calculates the excess masses associated with
                        the nodes of a level set tree
exmap                   Calculates a scale of excess mass profiles
hgrid                   Returns a grid of smoothing parameter values
leafsfirst              Calculates a level set tree, shape tree, or
                        tail tree
locofmax                Calculates the location of the maximum of a
                        function
lstseq.kern             Calculates a scale of kernel estimates
modecent                Returns locations of modes of a density
                        estimate
modegraph               Calculates a mode graph from a scale of
                        estimates
paracoor                Makes a parallel coordinate plot
pcf.func                Calculates a piecewise constant function for
                        some illustrative purposes
pcf.kern                Calculates a multivariate kernel estimate
plotbary                Makes a barycenter plot of a level set tree, a
                        location plot of a shape tree, or a tail tree
                        plot of a tail tree
plotbranchmap           Plots a branching map
plotexmap               Plots a scale of excess mass profiles
plotmodet               Plots a mode graph
plottree                Makes a tree plot of a level set tree, of a
                        shape tree, or of a tail tree
plotvolu                Makes a volume plot of a level set tree, a
                        shape plot of a shape tree, or a tail frequency
                        plot of a tail tree
plotvolu2d              Makes a perspective plot of a 2D volume
                        function or a 2D probability content function
profgene                Calculates the level set tree of a
                        rectangularwise constant function
profhist                Calculates the level set tree of a histogram
profkern                Calculates a level set tree of a kernel
                        estimate
prunemodes              Prunes modes away from a level set tree or a
                        shape tree
scaletable              Implements the scale and shape table
shape2d                 Returns a 2D volume function or 2D probability
                        content function
sim.data                Generates data for illustrative purposes
slicing                 Returns a one- or two-dimensional slice of a
                        multivariate function
stseq                   Calculates a sequence of radius functions from
                        a sequence of level sets
treedisc                Prunes a level set tree or a tail tree
tree.segme              Returns the segmentation of the nodes of 
                        a visualization tree
}

}
\author{
Jussi Klemelä <klemela@rumms.uni-mannheim.de>

Maintainer: Jussi Klemelä <klemela@rumms.uni-mannheim.de>

}

\references{ http://www.rni.helsinki.fi/~jsk/denpro }

\keyword{ multivariate }

\seealso{ 
\code{\link[delt]{delt}} 
}

\examples{
# level set tree

dendat<-sim.data(n=200,type="mulmod")   # data
pcf<-pcf.kern(dendat,h=1,N=c(32,32))    # kernel estimate

lst<-leafsfirst(pcf)             # level set tree
td<-treedisc(lst,pcf,ngrid=60)   # pruned level set tree

plotvolu(td)                     # volume plot

plotbary(td)                     # barycenter plot 

# shape tree

dendat<-sim.data(n=200,type="cross")    # data
pcf<-pcf.kern(dendat,h=1,N=c(32,32))    # kernel estimate

st<-leafsfirst(pcf,propo=0.01)   # shape tree, 1 per cent level set
tdst<-treedisc(st,pcf,ngrid=60)  # pruned shape tree

plotvolu(tdst)          # radius plot

plotbary(tdst)          # location plot

# tail tree

dendat<-sim.data(n=200,type="cross")    # data

tt<-leafsfirst(dendat=dendat,rho=0.65) # tail tree

plotbary(tt)     # tail tree plot

plotvolu(tt)     # tail frequency plot

}
