ccenuniapp<-function(uni){
#Calculates 1st moment of the union of rectangles "uni"
#
#uni is k*(2*d)-matrix, union of k rectangles in d-dimensional space.
#Returns a d-vector: 1st moment of uni: \int_uni x dx.
#
#Calls: centersit
#
k<-dim(uni)[1]
d<-dim(uni)[2]/2
depth<-3
valitulos<-matrix(0,d,depth)
if (k==1){             #if uni is only one rectangle
 tul<-centersit(uni)
 reslow=tul
 resupp=tul
}
else{
 endind<-seq(1:k)
 valitulos[,1]<-centersit(uni)
 cur<-uni
 for (taso in 2:min(depth,k)){
   curcur<-intersec(taso,endind,cur,uni)
   if (!is.na(curcur)){  #there is some intersection
     cur<-curcur$set
     endind<-curcur$endind  
     valitulos[,taso]<-centersit(cur) 
   }
 }
 reslow<-valitulos[,1]-valitulos[,2]
 resupp<-valitulos[,1]-valitulos[,2]+valitulos[,3]
}
return(list(low=reslow,hig=resupp))
}


