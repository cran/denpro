dist.func<-function(dendat,xepsi=0,yepsi=0,col="black",type="distr",
log="y",cex.axis=1)
{
n<-length(dendat)

if (type=="distr"){

plot(x="",y="",
xlim=c(min(dendat)-xepsi,max(dendat)+xepsi), ylim=c(0-yepsi,1+yepsi),
xlab="",ylab="",cex.axis=cex.axis)
ycur<-0
ordi<-order(dendat)
dendatord<-dendat[ordi]
for(i in 1:(n-1)){
    segments(dendatord[i],ycur,dendatord[i],ycur+1/n,col=col)
    segments(dendatord[i],ycur+1/n,dendatord[i+1],ycur+1/n,col=col)
    ycur<-ycur+1/n
}
segments(dendatord[n],ycur,dendatord[n],ycur+1/n,col=col)
segments(dendatord[n],ycur+1/n,max(dendat)+xepsi,ycur+1/n,col=col)

}
else if ((type=="right.tail") || (type=="left.tail")){

  if (type=="right.tail"){
         redu.ind<-(dendat>0) 
         dendat.redu<-dendat[redu.ind]
  }
  else{
         redu.ind<-(dendat<0)
         dendat.redu<--dendat[redu.ind]
  }
  ordi<-order(dendat.redu)
  dendat.ord<-dendat.redu[ordi]
  nredu<-length(dendat.redu)
  level<-seq(nredu,1)
  plot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis)

  #ordi<-order(dendat)
  #dendat.ord<-dendat[ordi]
  #medi.ind<-floor(n/2)
  #dendat.redu<-dendat.ord[medi.ind:n]
  #nredu<-length(dendat.redu)
  #level<-seq(nredu,1)
  #plot(dendat.redu,level,log="y",xlab="",ylab="")

}

}

