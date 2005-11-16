drawkern2d<-function(dendat,h,N,kernel="epane",trunc=3,threshold=0.0000001)
{
#
#source("~/kerle/profkernCRC.R")
#dyn.load("/home/jsk/kerle/kerCeva")
#dyn.load("/home/jsk/kerle/kerleCversio")
#pk2<-profkernCRC(dendat,h,N,Q)
#
#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(8,8)
#Q<-3

if (kernel=="gauss") h<-h*trunc   #trunc<-3

n<-dim(dendat)[1]
d<-length(N)
hnum<-length(h)
mnn<-maxnodenum(dendat,h,N,n,d)
extMaxnode<-mnn$maxnode
extMaxvals<-mnn$maxpositive
{
if (hnum>1){
 inh<-matrix(0,hnum+1,1)
 inh[2:(hnum+1)]<-h
}
else{
 inh<-h
}
}
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N

if (kernel=="epane") kertype<-1
else kertype<-2  # gaussian

kg<-.C("kergridC",
               as.integer(extMaxnode),
               as.integer(extMaxvals),
               as.double(dendat),
               as.double(inh),
               as.integer(inN),
               as.integer(n),
               as.integer(hnum),
               as.integer(d),
               as.integer(kertype),
               as.double(trunc), 
               as.double(threshold), 
               ioleft = integer(extMaxnode+1),
               ioright = integer(extMaxnode+1),
               ioparent = integer(extMaxnode+1),
               infopointer = integer(extMaxnode+1),
               iolow = integer(extMaxnode+1),
               ioupp = integer(extMaxnode+1),
               value = double(hnum*extMaxvals),
               index = integer(d*extMaxvals),
               nodefinder = integer(extMaxvals),
               numpositive = integer(1),
               numnode = integer(1),
PACKAGE="denpro")
#left<-kg$ioleft[2:(kg$numnode+1)]
#right<-kg$ioright[2:(kg$numnode+1)]
#parent<-kg$ioparent[2:(kg$numnode+1)]
#infopointer<-kg$infopointer[2:(kg$numnode+1)]
#iolow<-kg$iolow[2:(kg$numnode+1)]
#ioupp<-kg$ioupp[2:(kg$numnode+1)]

value<-kg$value[2:(kg$numpositive+1)]
#nodefinder<-kg$nodefinder[2:(kg$numpositive+1)]
vecindex<-kg$index[2:(d*kg$numpositive+1)]
index<-matrix(0,kg$numpositive,d)
for (i in 1:kg$numpositive){
  for (j in 1:d){
     index[i,j]<-vecindex[(i-1)*d+j]
  }
}

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

indenum<-dim(index)[1]

i<-1
while (i<=indenum){

   inde<-index[i,]
   point<-minim-h[1]+delta*inde
    
   z[1+inde[1],1+inde[2]]<-value[i]
   x[1+inde[1]]<-point[1]
   y[1+inde[2]]<-point[2]

   i<-i+1
}
x[1]<-minim[1]-h[1]
x[N[1]+2]<-minim[1]-h[1]+delta[1]*N[1]+delta[1]
y[1]<-minim[2]-h[1]
y[N[2]+2]<-minim[2]-h[1]+delta[2]*N[2]+delta[2]

return(x=x,y=y,z=z)
}









