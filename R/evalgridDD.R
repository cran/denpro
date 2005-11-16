evalgridDD<-function(M,sig,p,N,mul=3,theta=NULL,support=NULL)
# theta in 0:2*pi
{
d<-length(N)
mixnum<-length(p)
recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
  support<-matrix(0,2*d,1)
  for (i in 1:d){
      support[2*i-1]<-min(M[,i]-mul*sig[,i])
      support[2*i]<-max(M[,i]+mul*sig[,i])
  }
}

lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
    inde<-digit(i,N)
    if ((inde[1]==0) && (inde[2]==N[2])) inde<-c(0,0)
    point<-lowsuppo+step*inde-step/2
    if (!is.null(theta)){
         rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
         point<-rotmat%*%point
    }
    val<-0
    for (mi in 1:mixnum){
        evapoint<-(point-M[mi,])/sig[mi,]
        val<-val+p[mi]*evanor(evapoint)/prod(sig[mi,])
    }
    if (val>0){
       numpositive<-numpositive+1
       value[numpositive]<-val
       index[numpositive,]<-inde
    }
}

#index[2:recnum,]<-index[1:(recnum-1),]
#index<-index+1
#index[1,]<-c(1,1)

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index
high<-index+1

return(list(
value=value,index=index,
down=down,high=high,  #step=delta,
support=support,N=N))
}                              
