dfind<-function(dendat,ks){
#Finds for each observation max(ks) closest observations and
#distance to k:th closest observation for each ks[k].
#
#dendat is n*d-matrix
#ks levnum-vector of integers >=1
#
#Returns list(close,radius,dismat)
#close is n*max(ks)-matrix
#radius is n*levnum-matrix
#
n<-length(dendat[,1])
levnum<-length(ks)
radius<-matrix(0,n,levnum)
km<-max(ks)
#
#distan<-matrix(NA,n,n)
#for (i in 1:(n-1)){
#  for (j in (i+1):n){
#    eta<-etais(dendat[i,],dendat[j,])
#    distan[i,j]<-eta
#    distan[j,i]<-eta
#  }
#}
library(mva)
distan<-dist(dendat,method = "euclidean")
#
if (km==1){
  close<-matrix(0,n,1)
  for (i in 1:n){
      disrow<-as.matrix(distan)[i,]
      disrow[i]<-NA 
      obsind<-omaind(disrow)
      close[i]<-obsind
      radius[i,]<-disrow[obsind]/2 #we need dist only to the k:th closest
  }
}
else{
    close<-matrix(0,n,km)
    for (i in 1:n){
      disrow<-as.matrix(distan)[i,]
      disrow[i]<-NA          #NA is infinity for "omaind" 
      for (j in 1:km){
         obsind<-omaind(disrow)
         close[i,j]<-obsind
         for (l in 1:levnum){
            if (j==ks[l]) radius[i,l]<-disrow[obsind]/2
         }            #we need dist only to the ks[l]:th closest
         disrow[obsind]<-NA
      }
    }
}
return(list(close=close,radius=radius,dismat=distan))
}










