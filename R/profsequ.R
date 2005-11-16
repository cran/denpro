profsequ<-function(M,sig,p,h,N,Q)
#cvol=TRUE,ccen=TRUE,cfre=FALSE,compoinfo=FALSE)
{
mixnum<-dim(M)[1]
d<-dim(M)[2]

mul<-3
suppo<-matrix(0,2*d,1)

hnum<-length(h)
hrun<-1
while (hrun<=hnum){
   hcur<-h[hrun]

   sigkonvo<-sig+matrix(hcur,mixnum,d)   

   eg<-evalgrid(M,sigkonvo,p,N)

   for (i in 1:d){
     suppo[2*i-1]<-min(M[,i]-mul*sigkonvo[,i])
     suppo[2*i]<-max(M[,i]+mul*sigkonvo[,i])
   }

   curtree<-proftrem(eg,N,Q,suppo)

   if (hrun==1){
      if (hnum==1){
          treelist<-curtree
      }
      else{
          treelist=list(curtree)
      }
   }
   else{
      treelist=c(treelist,list(curtree))
   }
   hrun<-hrun+1
}
return(treelist)
}











