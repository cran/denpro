lststills<-function(pvec,M,sig,p,N,Q){

eva<-evalgrid(M,sig,p,N)

mul<-3
suppo<-matrix(0,2*d,1)
for (i in 1:d){
   suppo[2*i-1]<-min(M[,i]-mul*sig[,i])
   suppo[2*i]<-max(M[,i]+mul*sig[,i])
}
minim<-matrix(0,d,1)   #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-suppo[2*i-1]
  maxim[i]<-suppo[2*i]
}
delta<-(maxim-minim)/(N+1)
hmax<-0

dg<-drawdyak(eva$value,eva$index,N,minim,hmax,delta)

lenp<-length(pvec)

for (l in 1:lenp){

  p<-pvec[l]

  ######################

  big<-max(dg$z)
  lenx<-dim(dg$z)[1]
  leny<-dim(dg$z)[2]
  cutti<-p*big
  z2<-dg$z

  for (i in 1:lenx){
    for (j in 1:leny){
       if (dg$z[i,j]>=cutti) z2[i,j]<-cutti
    }
  }
  dg2<-dg
  dg2$z<-z2

  if (l==1) dglist<-list(dg2)
  else dglist<-c(dglist,list(dg2))

  #####################
  pr2<-proftrem(eva,N,Q,suppo,katka=p*big)

  pr2$level<-pr2$invalue

  if (l==1) prlist<-list(pr2)
  else prlist<-c(prlist,list(pr2))

}

return(list(dglist=dglist,prlist=prlist,big=big))
}


