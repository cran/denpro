leafsfirst.visu<-function(tt,et,lev=NULL,refe=NULL,type="lst",
levmet="radius",ordmet="etaisrec",
lkmbound=NULL,radius=NULL,
orde="furthest",suppo=T,propor=NULL,lty=NULL,numbers=TRUE,
sigcol="lightblue")
{

if ((!is.null(lev)) || (!is.null(propor))){
    type<-"shape"
    if (is.null(refe)) refe<-locofmax(et)
    if (!is.null(propor)) lev<-propor*max(pcf$value)
}
if (is.null(refe)) refe<-locofmax(et)

pp<-plotprof(tt,plot=FALSE,data=TRUE)
vecs<-pp$vecs

d<-length(et$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(et$support[2*i]-et$support[2*i-1])/et$N[i]

# order the atoms for the level set with level "lev"

lenni<-length(et$value)
distat<-matrix(0,lenni,1)
infopointer<-matrix(0,lenni,1)

if (type=="lst"){
  lkm<-lenni
  distat<-et$value
  infopointer<-seq(1,lkm)
}
else{

lkm<-0
for (i in 1:lenni){
  if (et$value[i]>=lev){
     lkm<-lkm+1
     nod<-i  #nod<-et$nodefinder[i]
     if (ordmet=="etaisrec"){
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-et$support[2*jj-1]+step[jj]*et$down[nod,jj]
            recci[2*jj]<-et$support[2*jj-1]+step[jj]*et$high[nod,jj]
         }
         distat[lkm]<-etaisrec(refe,recci)
     }
     else{
         lowi<-matrix(0,d,1)
         uppi<-matrix(0,d,1)
         for (jj in 1:d){
            lowi[jj]<-et$support[2*jj-1]+step[jj]*et$down[nod,jj]
            uppi[jj]<-et$support[2*jj-1]+step[jj]*et$high[nod,jj]
         }
         baryc<-lowi+(uppi-lowi)/2
         distat[lkm]<-etais(baryc,refe)  #etais(baryc[lk m,],baryind)
     }
     infopointer[lkm]<-i
  }
}

}  #else

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   #pointe->et$value,et$nodefinder

ord<-order(distat)
infopointer<-infopointer[ord]

if (suppo){
  xmin<-et$support[1]
  xmax<-et$support[2]
  ymin<-et$support[3]
  ymax<-et$support[4]
}
else{
  xmin<-tt$refe[1]-tt$maxdis  #et$support[1]
  xmax<-tt$refe[1]+tt$maxdis  #et$support[2]
  ymin<-tt$refe[1]-tt$maxdis  #et$support[3]
  ymax<-tt$refe[2]+tt$maxdis  #et$support[4]
}

plot(x=refe[1],y=refe[2],xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
pch=20) #,col="red")

i<-1
while (i<=lkm){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]   #ip<-et$nodefinder[infopointer[node]]

     x1<-et$support[1]+step[1]*et$down[ip,1]
     x2<-et$support[1]+step[1]*et$high[ip,1] 
     y1<-et$support[3]+step[2]*et$down[ip,2]
     y2<-et$support[3]+step[2]*et$high[ip,2] 
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="gray",lty=lty)

     i<-i+1
}

if (!is.null(lkmbound)){
  i<-1
  while (i<=lkmbound){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-et$nodefinder[infopointer[node]]

     x1<-et$support[1]+step[1]*et$down[ip,1]
     x2<-et$support[1]+step[1]*et$high[ip,1] 
     y1<-et$support[3]+step[2]*et$down[ip,2]
     y2<-et$support[3]+step[2]*et$high[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=sigcol,lty=lty)
     #points(x=refe[1],y=refe[2],pch=20,col="red")
     if (numbers) text(x=x1+(x2-x1)/2,y=y1+(y2-y1)/2,paste(i))

     i<-i+1
  }
}
else{
  i<-1
  radu<-tt$level[lkm]  #tt$madxdis
  while (radu>=radius){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-et$nodefinder[infopointer[node]]

     x1<-et$support[1]+step[1]*et$down[ip,1]
     x2<-et$support[1]+step[1]*et$high[ip,1] 
     y1<-et$support[3]+step[2]*et$down[ip,2]
     y2<-et$support[3]+step[2]*et$high[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="blue",lty=lty)
     points(x=refe[1],y=refe[2],pch=20,col="red")

     i<-i+1
     radu<-tt$level[node]
  }
}

}

