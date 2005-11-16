plot.kernscale<-function(scale,maxy0=0)
{
   hnum<-length(scale$hseq)
   for (i in 1:hnum){
     pk<-lk$pcfseq[[i]]
     dp<-draw.pcf(pk)
     if (i==1){ 
           minx<-min(dp$x) 
           miny<-min(dp$y)
           maxx<-max(dp$x) 
           maxy<-max(dp$y)  
     }
     else{ 
          minx<-min(minx,min(dp$x))
          miny<-min(miny,min(dp$y))
          maxx<-max(maxx,max(dp$x))
          maxy<-max(maxy,max(dp$y))
     }
   }
   maxy<-max(maxy,maxy0)
   plot(x="",y="",xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="",ylab="")
   for (i in 1:hnum){
     pk<-lk$pcfseq[[i]]
     dp<-draw.pcf(pk)
     matpoints(dp$x,dp$y,type="l")
   }
}


