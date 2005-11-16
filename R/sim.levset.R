sim.levset<-function(n,seed)
{
d<-2
mixnum<-11
D<-1.8
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)      #c(0,0)
M[2,]<-c(D,0)      #c(D1,0)
M[3,]<-c(2*D,0)
#M[4,]<-c(2.5*D,0)
M[4,]<-c(0,D)
M[5,]<-c(0,2*D)
M[6,]<-c(0,3*D)
M[7,]<-c(0,-D)
M[8,]<-c(0,-2*D)
M[9,]<-c(0,-3*D)
M[10,]<-c(-1.5,3.9*D)
M[11,]<-c(1.5,3.9*D)
sig<-matrix(1,mixnum,d)
sig[10,1]<-0.7
sig[11,1]<-0.7
p<-matrix(1,mixnum,1)
p[6]<-0.6
p[10]<-0.3
p[11]<-0.3
p<-p/sum(p)

dendat<-simmix(n,M,sig,p,seed=seed)

return(dendat)
}
