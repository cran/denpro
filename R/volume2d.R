volume2d<-function(radseq,gnum=500,type="radius",type2="slice",gnum2=1000,xmax=10)
{
# type "proba"    type2 "boundary"
lkm<-length(radseq$level)

if (type2=="slice"){

if (type=="radius") x<-radseq$level else x<-matrix(0,lkm,1)

td<-radseq$radseq[[1]]
if (type=="proba") td$volume<-td$proba
xy<-lst2xy(td,gnum=gnum)
ylen<-length(xy$x)
ystep<-1/(ylen-1)
y<-seq(0,1,ystep)   #matrix(0,xlen,1)
z<-matrix(0,length(x),length(y))

for (i in 1:lkm){
   td<-radseq$radseq[[i]]
   if (type=="proba"){ 
        tdvolume<-td$volume
        td$volume<-td$proba
        indi<-lkm-i+1
   }
   else indi<-i
   xy<-lst2xy(td,gnum=gnum)   #ma<-matchxy(xy$x,xy$y,y)

   if (type=="proba") x[indi]<-max(tdvolume) #[1]  #root=1

   ## normalize
   volu<-xy$x[length(xy$x)]-xy$x[1]
   int<-0
   step<-xy$x[2]-xy$x[1]
   for (j in 1:length(xy$x)){
       int<-int+step*xy$y[j]
   }
   b<-volu^2/int
   ynew<-b*xy$y
   ## end normalize
   z[indi,]<-ynew   
}
}
else{ #type2=="boundary"

if (is.null(xmax)){
    td<-radseq$radseq[[1]]
    if (type=="proba") td$volume<-td$proba
    xmax<-max(td$volume)
}

ymax<-xmax
step<-2*xmax/(gnum-1)
x<-seq(-xmax,xmax,step)
y<-x
z<-matrix(0,length(x),length(y))

for (i in 1:lkm){
  td<-radseq$radseq[[i]]
  if (type=="proba") td$volume<-td$proba
  xy<-lst2xy(td,gnum=gnum2,type=type)  

  ## normalize
  volu<-xy$x[length(xy$x)]-xy$x[1]
  int<-0
  step<-xy$x[2]-xy$x[1]
  for (j in 1:length(xy$x)){
       int<-int+step*xy$y[j]
  }
  b<-volu^2/int
  ynew<-b*xy$y
  ## end normalize

  for (j in 1:length(x)){
      for (k in 1:length(y)){
          len<-sqrt(x[j]^2+y[k]^2)
          xn<-x[j]/len
          yn<-y[k]/len
          th2<-atan(xn/yn)
          if (yn<0) th2<-atan(xn/yn)+pi else if (xn<0) th2<-atan(xn/yn)+2*pi
          propo<-th2/(2*pi) 
          dirind<-max(1,round( propo*length(xy$x) ))
          rho<-ynew[dirind]
          if (len<=rho) z[j,k]<-radseq$level[i]
      }
  }
}

}

return(list(x=x,y=y,z=z,type=type,type2=type2))
}


