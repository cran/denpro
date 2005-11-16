drawmix<-function(M,sig,p,plkm){
#Makes data for drawing a mixture of normal densities
#
#M is mixnum*2-matrix, means
#sig is mixnum*2-matrix, variances
#p is mixnum-vector
#plkm is large integer, 30
#
minmx<-min(M[,1])
maxmx<-max(M[,1])
minmy<-min(M[,2])
maxmy<-max(M[,2])
maxsx<-max(sig[,1])
maxsy<-max(sig[,2])
#
xnum<-plkm
ynum<-plkm
#
xbeg<-minmx-4*maxsx
xend<-maxmx+4*maxsx
ybeg<-minmy-4*maxsy
yend<-maxmy+4*maxsy
xstep<-(xend-xbeg)/(xnum-1)
ystep<-(yend-ybeg)/(ynum-1)
xseq<-seq(xbeg,xend,xstep)
yseq<-seq(ybeg,yend,ystep)
#
mixnum<-length(p)
z<-matrix(0,plkm,plkm)
for (i in 1:mixnum){
  z<-z+p[i]*drawnor(M[i,],sig[i,],xseq,yseq)
}
return(list(x=xseq,y=yseq,z=z))
}
