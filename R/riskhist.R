riskhist<-function(frekv,recs,n){
#Estimates expected integrated squared error minus \int f^2
#
#frekv is binlm-vector, frequency for every bin
#recs is binlkm*(2*d)-matrix, bins
#
#Returns a real number
#
#Estimate is \sum_{i=1}^binlkm (1-ni*(1+1/n))*ni/(n^2*tili)
#
if (dim(t(recs))[1]==1) binlkm<-1 else binlkm<-length(recs[,1]) 
if (binlkm==1){
  ans<--1/(massone(recs))
}
else{ 
    sum<-0
    for (i in 1:binlkm){
        tili<-massone(recs[i,])
        ni<-frekv[i]
        sum<-sum+(1-ni*(1+1/n))*ni/(n^2*tili)
    }
    ans<-sum
}
return(ans)
}
