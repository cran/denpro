sim.data<-function(n=NULL,seed=1,N=NULL,type="mulmod",
M=NULL,sig=NULL,p=NULL,d=NULL,
cova=NULL,marginal="student",t=NULL,df=NULL,distr=FALSE, noisedim=1,
sig1=0.5,sig2=1.5,diff=0.1
)
{
if (type=="mixt") return( simmix(n,M,sig,p,seed,d) )

if (type=="mulmod") return( sim.mulmod(n=n,seed=seed,N=N) )

if (type=="fox") return( sim.fox(n=n,seed=seed,N=N) )

if (type=="tetra3d") return( sim.tetra3d(n=n,seed=seed,N=N) )

if (type=="penta4d") return( sim.penta4d(n=n,seed=seed,N=N) )

if (type=="cross") return( sim.cross(n=n,seed=seed,N=N,sig1=sig1,sig2=sig2) )

if (type=="1d2modal") return( sim.1d2modal(n=n,seed=seed,N=N,distr=distr) )

if (type=="claw") return( sim.claw(n=n,seed=seed,N=N) )

if (type=="fssk") return( sim.fssk(n=n,noisedim=noisedim,seed=seed) )

if (type=="nested") return( sim.nested(n=n,seed=seed,N=N) )

if (type=="mulmodII") return( sim.mulmodII(n=n,seed=seed,N=N) )

if (type=="gauss"){
   eig<-eigen(cova,symmetric=TRUE)
   sigsqm<-eig$vectors%*%diag(eig$values^{1/2})  
   set.seed(seed)
   symmedata<-matrix(rnorm(2*n),n,2)
   dendat<-t(sigsqm%*%t(symmedata))
   dendat<-pnorm(dendat)

   if (marginal=="gauss") dendat<-qnorm(dendat)
   if (marginal=="student") dendat<-qt(dendat, df=t)
   return(dendat)
}

if (type=="student"){
   eig<-eigen(cova,symmetric=TRUE)
   sigsqm<-eig$vectors%*%diag(eig$values^{1/2})  
   set.seed(seed)
   symmedata<-matrix(rt(2*n,df=df),n,2)
   dendat<-t(sigsqm%*%t(symmedata))
   dendat<-pt(dendat,df=df)

   if (marginal=="gauss") dendat<-qnorm(dendat)
   if (marginal=="student") dendat<-qt(dendat, df=t)
   return(dendat)
}

if (type=="gumbel"){
  link<-function(y,g){ return ( (-log(y))^g ) }
  linkinv<-function(y,g){ return ( exp(-y^(1/g)) ) }
  der1<-function(y,g){ return ( -g*(-log(y))^(g-1)/y ) }
  der1inv<-function(y,g){ return ( y ) }
}

if (type=="diff1d"){
   xala<--0
   xyla<-1
   support<-c(xala,xyla)
   d<-1
   M<-c(0.5-diff,0.5+diff)
   mixnum<-length(M)
   sig<-matrix(sig1,mixnum,d)
   p<-matrix(1,mixnum,1)
   p<-p/sum(p)
   pcf<-pcf.func("mixt",N=N,sig=sig,M=M,p=p,support=support,distr=distr)
   return(pcf)
}

}

