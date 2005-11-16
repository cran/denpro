drawkern<-function(value,index,N,
dendat=NULL,h=NULL)
{

d<-length(N)

if (d==2){

x<-matrix(0,N[1]+2,1)
y<-matrix(0,N[2]+2,1)
z<-matrix(0,N[1]+2,N[2]+2)

minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-min(dendat[,i])
  maxim[i]<-max(dendat[,i])
}
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)
hmax<-h[1]

indenum<-dim(index)[1]

i<-1
while (i<=indenum){

   inde<-index[i,]
   #point<-minim-h[1]+delta*inde

   z[1+inde[1],1+inde[2]]<-value[i]
   #x[1+inde[1]]<-point[1]
   #y[1+inde[2]]<-point[2]

   i<-i+1
}

i<-1
while (i<=N[1]){
   x[1+i]<-minim[1]-hmax+delta[1]*i
   i<-i+1
}

i<-1
while (i<=N[2]){
   y[1+i]<-minim[2]-hmax+delta[2]*i
   i<-i+1
}

x[1]<-minim[1]-h[1]
x[N[1]+2]<-minim[1]-h[1]+delta[1]*N[1]+delta[1]
y[1]<-minim[2]-h[1]
y[N[2]+2]<-minim[2]-h[1]+delta[2]*N[2]+delta[2]

return(x=x,y=y,z=z)

}

else{    #d=1


x<-matrix(0,N+2,1)
y<-matrix(0,N+2,1)

minim<-min(dendat)
maxim<-max(dendat)
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)

indenum<-dim(index)[1]

i<-1
while (i<=indenum){

   inde<-index[i]
   point<-minim-hmax+delta*inde
    
   y[1+inde]<-value[i]
   x[1+inde]<-point

   i<-i+1
}
x[1]<-minim-hmax
x[N+2]<-minim-hmax+delta*N+delta

return(x=x,y=y)

}

}











