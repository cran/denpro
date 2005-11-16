qq.plot<-function(dendat=NULL,compa="gauss",basis="gauss",
mean=0,sig=1,df=1,
gnum=1000,d=1,R=3,qqtype="1d")
{
if (qqtype=="1d"){
   n<-dim(dendat)[1]
   p<-(seq(1:n)-1/2)/n
   if (compa=="gauss") x<-qnorm(p,mean=mean,sd=sig)
   if (compa=="student") x<-qt(p,df=df)
   if (compa=="unif") x<-qunif(p)
   if (compa=="exp") x<-qexp(p)
   y<-dendat[order(dendat)]
   tyyppi<-"p"
   xlab<-"compared quantiles"
   ylab<-"empirical quantiles"
}
if (qqtype=="p2v"){
     rp<-tailfunc(R,d,type=compa,gnum=gnum,sig=sig,nu=df)
     x<-rp$volu
     rp2<-tailfunc(R,d,type=basis,gnum=gnum,sig=sig,nu=df)
     y<-rp2$volu
     tyyppi="l"
     ylab<-"empirical"
     xlab<-"model"
}

plot(x,y,type=tyyppi,ylab=ylab,xlab=xlab)
maxxy<-max(max(x),max(y))
minxy<-min(min(x),min(y))
segments(minxy,minxy,maxxy,maxxy)
}





