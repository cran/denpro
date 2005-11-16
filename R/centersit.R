centersit<-function(rec){
#Calculates a sum of 1st moments of a set of rectangles
#
#rec is k*(2*d)-matrix, represents k rectangles in d-space
#We may assume that rec[j,2i]>rec[j2i-1]
#
#Returns a d-vector, elements are \int_R x_jdx, where R
#is k:th row of rec and j=1,...,d.
#
if (dim(t(rec))[1]==1) k<-1 else k<-length(rec[,1])  #rows of rec
if (k==1){
 d<-length(rec)/2
 tulos<-matrix(0,d,1)
 j<-1
 while (j<=d){
    apurec<-rec      #apurec such that is volume is equal to
    apurec[2*j-1]<-0 #volume of d-1 dimensional rectangle where
    apurec[2*j]<-1   #we have removed j:th dimension
    vajmas<-massat(apurec) 
    tulos[j]<-vajmas*(rec[2*j]^2-rec[2*j-1]^2)/2
    j<-j+1
 }
}
else{
 d<-length(rec[1,])/2
 tulos<-matrix(0,d,1)
 mattulos<-matrix(0,k,d)
 for (i in 1:k){
   j<-1
   while (j<=d){
     apurec<-rec[i,]  #apurec such that is volume is equal to
     apurec[2*j-1]<-0 #volume of d-1 dimensional rectangle where
     apurec[2*j]<-1   #we have removed j:th dimension
     vajmas<-massat(apurec) 
     mattulos[i,j]<-vajmas*(rec[i,2*j]^2-rec[i,2*j-1]^2)/2
     j<-j+1
   }
   tulos<-tulos+mattulos[i,]
 }
}
return(tulos)
}


