modeorder<-function(tree){
#
if (modlkm>1){
  modesInOrder<-matrix(0,modlkm,1)  #pointers to modes, in right order
  if (ordercrit==0){  
     if (is.null(crit)){
        d<-dim(center)[1]     #dimension is the number of rows of center
        crit<-matrix(0,d,1)
     }
     #
     #order so that on the leftmost is the closest to crit
     # after that 2 possibilities: 
     #1. second is second closest to origo
     #2. second is closest to first  
     #we choose 2. 
     #
     distanceToOrigo<-matrix(NA,modlkm,1)    #NA is infty 
     for (i in 1:modlkm){
         cur<-modloc[i]
         curre<-center[,cur]
         distanceToOrigo[i]<-etais(curre,crit)
     }
     modesInOrder<-omaord2(modloc,distanceToOrigo)  
  }
  #else 
  #order so that on the left is mode with the smallest "ordercrit"-coordinate
}
return(modesInOrder)
}

