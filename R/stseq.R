stseq<-function(N,lnum,
func=NULL,dendat=NULL,
h=NULL,Q=NULL,kernel="epane",hw=NULL,
sig=rep(1,length(N)),support=NULL,theta=NULL,
M=NULL,p=NULL,mul=3,
t=1,marginal="unif",r=0,
mu=NULL,xi=NULL,Omega=NULL,alpha=NULL,df=NULL,g=1
)
{
#lnum<-length(lseq)
level<-matrix(0,lnum,1)
volume<-matrix(0,lnum,1)
if (!is.null(dendat)) pcf<-pcf.kern(dendat,h,N,kernel=kernel,hw=hw)
else
pcf<-eval.func(func,N,
sig=sig,support=support,theta=theta,
M=M,p=p,mul=mul,
t=t,marginal=marginal,r=r,
mu=mu,xi=xi,Omega=Omega,alpha=alpha,df=df,g=g
)

maksi<-max(pcf$value)
for (i in 1:lnum){   
      lev<-maksi*i/(lnum+1) 
      level[i]<-lev 
      refe<-locofmax(pcf)
      st<-leafsfirst(pcf,lev=lev,refe=refe)
      volume[i]<-max(st$volume)
      if (i==1){
           if (lnum==1){ 
               istseq<-st
           }
           else{
               stseq<-list(st)
           }
      }
      else{
          stseq<-c(stseq,list(st))
      }
}
return(list(radseq=stseq,level=level,volume=volume))
}



