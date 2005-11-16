drawdyak<-function(value,index,N,
minim=NULL,hmax=NULL,delta=NULL){
#
#ke<-kergrid(dendat,h,N) 
#tul<-drawdyak(ke$value,ke$index,N)
#persp(tul$x,tul$y,tul$z,,theta=130,phi=20) 
#plot(tul$x,tul$z)
#
d<-length(N)
#
if (d==2){
  xdim<-N[1]
  ydim<-N[2]
  x<-seq(1,xdim+2)
  y<-seq(1,ydim+2)
  z<-matrix(0,xdim+2,ydim+2)
  posnum<-length(value)

  if (!is.null(minim) && !is.null(hmax) && !is.null(delta)){
     for  (ip in 1:(xdim+2)){
        x[ip]<-minim[1]+(ip-1)*delta[1]
     }

     for  (ip in 1:(ydim+2)){
        y[ip]<-minim[2]+(ip-1)*delta[2]
     }
  }
 
  for (i in 1:posnum){
     xcoor<-index[i,1]
     ycoor<-index[i,2]
     z[xcoor+1,ycoor+1]<-value[i]   
  }

}
else{  #d=1
  xdim<-N
  x<-seq(1,xdim+2)
  z<-matrix(0,xdim+2,1)
  posnum<-length(value)
  for (i in 1:posnum){

    if (!is.null(minim) && !is.null(hmax) && !is.null(delta)){
        inde<-index[i]
        point<-minim-hmax+delta*inde
        x[1]<-minim-hmax
        x[xdim+2]<-minim-hmax+delta*N+delta
    }
    else point<-index[i]

    xcoor<-index[i]
    x[xcoor+1]<-point
    z[xcoor+1]<-value[i]   
  }
  y<-NULL
}

return(list(x=x,y=y,z=z))
}






