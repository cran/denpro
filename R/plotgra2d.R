plotgra2d<-function(pre,dendat,d1,d2,zdim,xlab,ylab){
#
#seplsets,newradius,sepval
#
#seplsets is componum*n-matrix
#newradius is n*componum-matrix
#sepval is componum-vector
#
seplsets<-pre$seplsets
newradius<-pre$newradius
sepval<-pre$sepval
#
n<-dim(dendat)[1]
dendat2d<-matrix(0,n,2)
dendat2d[,1]<-dendat[,d1]
dendat2d[,2]<-dendat[,d2]
#
z<-matrix(0,zdim,zdim)
componum<-dim(newradius)[2]
#
minx<-min(dendat2d[,1])
maxx<-max(dendat2d[,1])
miny<-min(dendat2d[,2])
maxy<-max(dendat2d[,2])
maxrad<-max(newradius)
eta<-0.1
#
xbeg<-minx-maxrad-eta
xend<-maxx+maxrad+eta
ybeg<-miny-maxrad-eta
yend<-maxy+maxrad+eta
#
xdelta<-(xend-xbeg)/(zdim+1)
ydelta<-(yend-ybeg)/(zdim+1)
#
#[i,j] on xbeg+i*xdelta
#
for (i in 1:componum){
   for (j in 1:n){
      if (seplsets[i,j]==1){  
         obse<-dendat2d[j,]
         ini<-round((obse[1]-xbeg)/xdelta)
         inj<-round((obse[2]-ybeg)/ydelta)
         radcur<-newradius[j,i]
         xh<-floor(radcur/xdelta)
         yh<-floor(radcur/ydelta)
         for (k in (ini-xh):(ini+xh)){
            for (l in (inj-yh):(inj+yh)){
               xcur<-xbeg+k*xdelta
               ycur<-ybeg+l*ydelta
               if (((xcur-obse[1])^2+(ycur-obse[2])^2)<=radcur^2){
                    z[k,l]<-sepval[i]
               }
            }
         }
      }
   }
}
colonum<-componum+1
#colo<-rainbow(colonum,s=1, v=1, start=0, end=max(1,colonum - 1)/n, gamma=1)
#colo<-heat.colors(colonum)
colo<-terrain.colors(colonum)
#colo<-topo.colors(colonum)
#colo<-cm.colors(colonum)
colo[1]<-"white"
#
xlen<-zdim
x<-matrix(0,xlen,1)
ylen<-zdim
y<-matrix(0,ylen,1)
x[1]<-xbeg
for (i in 2:xlen){
   x[i]<-x[i-1]+xdelta
}
y[1]<-ybeg
for (i in 2:ylen){
   y[i]<-y[i-1]+ydelta
}                           
#
image(x,y,z,col=colo,xlab=xlab,ylab=ylab)    #add=T
points(dendat2d,pch=20)            
}







