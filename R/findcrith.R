findcrith<-function(mt,snum)
{
crith<-matrix(0,snum,1)

lenni<-length(mt$ycoor)
mlabel<-mt$mlabel
low<-matrix(0,snum,1)
upp<-matrix(0,snum,1)
low[1]<-1
glob<-2
while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
}
upp[1]<-glob-1
# now glob is at the start of new block
i<-2
while (i<=snum){
   low[i]<-glob
   glob<-glob+1
   while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
   }
   upp[i]<-glob-1
   i<-i+1
}

prenum<-upp[1]-low[1]+1     #number of modes
lkm<-0
i<-2
while (i<=snum){
  curnum<-upp[i]-low[i]+1
  if (curnum!=prenum){
        lkm<-lkm+1
        crith[lkm]<-i
  }
  i<-i+1
  prenum<-curnum
}
crith<-crith[1:lkm]

return(crith)
}
