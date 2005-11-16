ccenuni<-function(uni){
#Calculates 1st moment of the union of rectangles "uni"
#
#uni is k*(2*d)-matrix, union of k rectangles in d-dimensional space.
#Returns a d-vector: 1st moment of uni: \int_uni x dx.
#
#Calls: centersit
#
if (dim(t(uni))[1]==1) k<-1 else k<-length(uni[,1]) #rows of uni
if (k==1){             #if uni is only one rectangle
 tul<-centersit(uni)
 #d<-length(uni)/2
 #tul<-1
 #  for (j in 1:d){
 #    tul<-tul*(uni[2*j]-uni[2*j-1])
 #  }
}
else{
 endind<-seq(1:k)
 tul<-centersit(uni)
 cur<-uni
 for (taso in 2:k){
   curcur<-intersec(taso,endind,cur,uni)
   if (!is.na(curcur)){  #there is some intersection
     cur<-curcur$set
     endind<-curcur$endind  
     cent<-centersit(cur) 
     tul<-tul+(-1)^(taso-1)*cent
   }
 }
}
return(tul)
}


