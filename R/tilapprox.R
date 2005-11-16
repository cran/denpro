tilapprox<-function(runi){
#
lkm<-dim(runi)[1]   #number of rows
masses<-matrix(0,3,1)  #we calculate sum of masses, sum of masses of
                       #pairwise intersectiosn, and sum of three wise
#
masses[1]<-sum(massat(runi))    #kaiteiden massojen summa
#
if (lkm>=2){           #if more than 1 rectangle
 apu<-til1(runi)
 ind<-apu$ind
 curkosk<-apu$curkosk
 currecs<-apu$currecs
 parimat<-apu$parimat
 if (ind>0){ #jos oli parittaisia leikkauksia
    masses[2]<-sum(massat(currecs)) #parittaisten leikkausten massojen summa
    kosk<-3
    apu2<-til2(runi,curkosk,currecs,parimat,kosk)
    ind<-apu2$ind
    if (ind>0){
       currecs<-apu2$currecs
       curkosk<-apu2$curkosk
       masses[kosk]<-sum(massat(currecs))
    }
 }
}
reslow=masses[1]-masses[2]
resupp=masses[1]-masses[2]+masses[3]
return(list(reslow=reslow,resupp=resupp))
}
