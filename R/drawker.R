drawker<-function(dendat,h,hillkm){
#
#persp(tul$x,tul$y,tul$z,phi=20)   
#
sade<-2
xmax<-max(dendat[,1])+sade
xmin<-min(dendat[,1])-sade
ymax<-max(dendat[,2])+sade
ymin<-min(dendat[,2])-sade
xaskel<-(xmax-xmin)/(hillkm-1)
yaskel<-(ymax-ymin)/(hillkm-1)
x<-seq(xmin,xmax,xaskel) 
y<-seq(ymin,ymax,yaskel)
lx<-length(x)
ly<-length(y)
z<-matrix(0,lx,ly)
for (i in 1:lx){
   for (j in 1:ly){
     z[i,j]<-kerest2d(x[i],y[j],h,dendat) 
   }
}
return(list(x=x,y=y,z=z))
}

