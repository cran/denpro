fstbranch<-function(mt){
#Finds the part of the lev.sets which branches first
#
#mt is a multitree
#
roots<-mt$roots
child<-mt$child
sibling<-mt$sibling
#
firstbranch<-NA
#
prrootnum<-length(roots)
if (prrootnum>1){
   firstbranch<-roots[1]
}
else{
      cur<-roots[1]
      while (child[cur]>0){
         cur<-child[cur]
         if (sibling[cur]>0){       
            firstbranch<-cur 
         }
      }
} 
return(firstbranch)
}

