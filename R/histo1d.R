histo1d<-function(dendat,binlkm,ala=NULL,yla=NULL,
pic=TRUE,brush=NULL,brushcol=c("blue"),col=NULL,border=NULL,
xlab="",ylab="",cex.lab=1,cex.axis=1,data=FALSE)
{
if (is.null(ala)) ala<-min(dendat)
if (is.null(yla)) yla<-max(dendat)
step<-(yla-ala)/binlkm
frekv<-matrix(0,binlkm,1)
value<-matrix(0,binlkm,1)
if (!is.null(brush)){
   cnum<-max(brush)
   shade <-matrix(0,binlkm,cnum)
}
n<-length(dendat)
for (i in 1:n){
   hava<-dendat[i]
   ind<-min(binlkm,floor((hava-ala)/step)+1)
   frekv[ind]<-frekv[ind]+1
   if ((!is.null(brush)) && (brush[i]>0)) 
              shade[ind,brush[i]]<-shade[ind,brush[i]]+1
}
value<-frekv/(n*step)

if (pic){
   plot(x="",y="",xlab=xlab,ylab=ylab,xlim=c(ala,yla),ylim=c(0,max(value)),
   cex.lab=cex.lab,cex.axis=cex.axis)
   for (i in 1:binlkm){
          xala<-ala+(i-1)*step
          xyla<-xala+step
          y<-value[i]
 
          polygon(c(xala,xala,xyla,xyla),c(0,y,y,0),col=col,border=border)

          if (!is.null(brush)){
              y0<-0
              for (kk in 1:cnum){
                  y<-y0+shade[i,kk]
                  polygon(c(xala,xala,xyla,xyla),c(y0,y,y,y0),col=brushcol[kk])
                  y0<-y
              }
          }
      }
}
if (data){
     return(list(frekv=frekv,ala=ala,step=step,value=value))
}
}



