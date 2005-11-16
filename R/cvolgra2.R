cvolgra2<-function(seplsets,newradius,dendat,approx=TRUE){
#Laskee joukolle tasojoukon osia voluumin.
#
#"seplsets" is tasolkm*n-matriisi, 0:a ja 1:a
#"newradius" is n*tasolkm-matrix, newradiuses separately for each level 
#"dendat" is n*d-matrix
#
#Returns a tasolkm-vector: volume for each set.
#            
n<-length(dendat[,1])
radiushig<-newradius
radiuslow<-newradius/sqrt(2)
tasolkm<-length(seplsets[,1])     #seplsets:n rivien maara
volshig<-matrix(0,tasolkm,1)
volshigupp<-matrix(0,tasolkm,1)
volslow<-matrix(0,tasolkm,1)
volslowupp<-matrix(0,tasolkm,1)
for (a in 1:tasolkm){
  radsh<-radiushig[,a]
  runi<-rec2rec(seplsets[a,],radsh,dendat)
  if (approx){
     appu<-tilapprox(runi)
     alla<-appu$reslow
     ylla<-appu$resupp
     volshig[a]<-alla
     volshigupp[a]<-ylla
  }
  else volshig[a]<-til(runi)
}
for (a in 1:tasolkm){
   radsl<-radiuslow[,a]
   runi<-rec2rec(seplsets[a,],radsl,dendat)
   if (approx){
     appu<-tilapprox(runi)
     alla<-appu$reslow
     ylla<-appu$resupp
     volslow[a]<-alla
     volslowupp[a]<-ylla
   }
   else volslow[a]<-til(runi)
}
#vols<-(volshig+volslow)/2
#return(t(vols))
return(list(lowlow=volslow,lowhig=volslowupp,higlow=volshig,highig=volshigupp))
}             





