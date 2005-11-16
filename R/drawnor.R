drawnor<-function(mu,s,xseq,yseq){
#Makes data for drawing perspective plot from normal-density function
#
#mu is 2-vector, means
#s is 2-vector, variances
#plkm is large integer
#
#nor<-drawnor(mu,s,plkm)
#persp(nor$x,nor$y,nor$z) 
#
xnum<-length(xseq)
ynum<-length(yseq)
#
X<-xseq*matrix(1,xnum,xnum)
Y<-t(yseq*matrix(1,ynum,ynum))
z<-exp(-(X-mu[1])^2/(2*s[1]^2)-(Y-mu[2])^2/(2*s[2]^2))/(2*pi*s[1]*s[2])
#
return(z)
}

