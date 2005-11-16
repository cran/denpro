drawmix1d<-function(mu,sig,p,plkm){
#Makes data for drawing a mixture of normal densities
#
#mu is mixnum-matrix, means
#sig is mixnum-matrix, variances
#p is mixnum-vector
#plkm is large integer, 30
#
minmu<-min(mu)
maxmu<-max(mu)
minsig<-min(sig)
maxsig<-max(sig)
#
xbeg<-minmu-4*maxsig
xend<-maxmu+4*maxsig
xstep<-(xend-xbeg)/(plkm-1)
xseq<-seq(xbeg,xend,xstep)
#
mixnum<-length(p)
z<-matrix(0,plkm,plkm)
for (i in 1:mixnum){
  z<-z+p[i]*drawnor1d(mu[i],sig[i],xseq)
}
return(list(x=xseq,z=z))
}
