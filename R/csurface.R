csurface<-function(component,AtomlistNext,AtomlistAtom,
minim,h,delta,index,d,center){

componum<-length(component)
radius<-matrix(0,componum,2*d)

for (i in 1:componum){
   pointer<-component[i]
   barycenter<-center[,i]
   radius[i,1:d]<-barycenter
   radius[i,(d+1):(2*d)]<-barycenter
   while (pointer>0){
        atompointer<-AtomlistAtom[pointer]
        inde<-index[atompointer,]
        newcente<-minim-h+delta*inde
        newlow<-newcente-delta/2
        newupp<-newcente+delta/2
        for (k in 1:d){
           if (newlow[k]<=radius[i,k]) radius[i,k]<-newlow[k]
           if (newupp[k]>=radius[i,d+k]) radius[i,d+k]<-newupp[k]
        }
        pointer<-AtomlistNext[pointer]
   }
}
return(radius)
}
