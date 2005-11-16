drawnor1d<-function(mu,s,xseq){
#Makes data for drawing perspective plot from normal-density function
#
#mu is mean
#s is variance
#plkm is large integer
#
#z<-drawnor1d(mu,s,plkm)
#plot(xseq,z) 
#
xnum<-length(xseq)
#
X<-t(xseq*matrix(1,xnum,1))
z<-exp(-(X-mu)^2/(2*s))/sqrt(2*pi)
#
return(z)
}

