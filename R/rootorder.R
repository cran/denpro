rootorder<-function(roots,siborder){
#
rootnum<-length(roots)
rord<-matrix(0,rootnum,1)
for (i in 1:rootnum){
   r<-roots[i] 
   ind<-siborder[r]
   rord[ind]<-r
}
return(rord)
}
