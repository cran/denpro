makeker2<-function(eval,levnum,ks,close){
#Represents kernel estimate with level sets
#
#eval is n-vector, evaluation of a kernel estimate at datapoints
#levnum is the number of desired level sets
#ks is levnum-vector: numbers of nearest neighbours to be taken into account
#  in determining the radiuses of balls centered at datapoints
#  (level sets are unions of balls centered at datapoints)
#close is n*max(ks)-matrix, max(ks) closest observations for each observation
#radius is n*levnum-matrix, radius for each observation, in each level
#
#Returns list(values,lsets,max)
# values is levnum-vector
# lsets is levnum*n-matrix
# max>0, tarvitaan piirtoa varten
#
n<-length(eval)
tieto<-matrix(0,levnum,1)
lsets<-matrix(0,levnum,n)      
lsets2<-matrix(NA,levnum,n)
#
M<-max(eval)  #tih f:on maksimi estimoidaan laskemalla maksimi hav pist:ssa
#askel<-M/levnum       #tasot ovat 0,askel,...,M-askel
#tieto<-t(t(seq(alin,M-askel,askel)))  #pal.vakion f:on arvot tietoon
askel<-M/(levnum+1)       #tasot ovat askel,2*askel,...,M-askel
tieto<-t(t(seq(askel,M-askel,askel)))  #pal.vakion f:on arvot tietoon
#
curcol<-matrix(1,levnum,1)
for (i in 1:n){
  val<-eval[i]
  ind<-floor(val/askel)
  if (ind==0) ind<-1
  if (ind>levnum) ind<-levnum
  lsets[1:ind,i]<-1
  for (j in 1:ind){
    ccol<-curcol[j]
    lsets2[j,ccol]<-i
    curcol[j]<-ccol+1
  }
}
#kaydaan tasot lapi ja poistetaan ne tasojoukon
#havainnot, jotka eivat riittavasti erillaan komplementista  
notempty<-TRUE  #highest lsetsets may become empty
empind<-0
i<-1
while ((i<=levnum) && (notempty)){            
  j<-1
  while ((j<=n) && (!is.na(lsets2[i,j]))){  
    curobs<-lsets2[i,j]
    for (l in 1:ks[i]){
       nearby<-close[curobs,l]
       if (lsets[i,nearby]==0)  lsets[i,curobs]<-0
           #if observation close to curobs is on the complement of
           #the i:th level-set, we remove this obs 
    }
    j<-j+1
  }
  if (sum(lsets[i,])==0){ 
        notempty<-FALSE
        empind<-empind+1
  }
  i<-i+1
}
notempty<-levnum-empind
tieto<-t(t(tieto[1:notempty]))  #save only nonempty
lsets<-lsets[1:notempty,]
#ks<-t(t(ks[1:notempty]))
return(list(levels=tieto,lsets=lsets,max=M))
}












