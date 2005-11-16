evalgrid1D<-function(M,sig,p,N,mul=3)
{
mixnum<-length(p)
value<-matrix(0,N,1)
index<-seq(1,N)

support<-matrix(0,2,1)
support[1]<-min(M-mul*sig)
support[2]<-max(M+mul*sig)
step<-(support[2]-support[1])/N

for (i in 1:N){
    point<-support[1]+step*index[i]-step/2
    val<-0
    for (j in 1:mixnum){
        val<-val+p[j]*evanor(point-M[j]/sig[j])/sig[j]
    }
    value[i]<-val
}

return(list(
value=value,index=index,
#down=down,high=high,step=delta,
support=support,N=N))
}                              
