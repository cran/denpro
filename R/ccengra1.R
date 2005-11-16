ccengra1<-function(levels,radius,dendat,maslow,mashig){
#Laskee joukolle tasojoukon osia painopisteen. 
#
#"levels" is tasolkm*n-matriisi, 0:a ja 1:a
#"radius" is tasolkm-vector
#"dendat" is n*d-matrix
#maslow, mashig are tasolkm-vectors: 
#masses of sets for radius=radiushig, radius=radiuslow correspondingly
#
#Returns a tasolkm*d-matrix: center of mass for each set.
#
#calls: ccenuni
#
n<-length(dendat[,1]) 
radiushig<-radius
radiuslow<-radius/sqrt(2)
tasolkm<-length(levels[,1])     #levels:n rivien maara
d<-length(dendat[1,])
volshig<-matrix(0,tasolkm,d)
volslow<-matrix(0,tasolkm,d)
for (a in 1:tasolkm){
  radsh<-radiushig[a]*matrix(1,n,1)    
  uni<-rec2rec(levels[a,],radsh,dendat)   
  volshig[a,]<-ccenuni(uni)/mashig[a]
}
for (a in 1:tasolkm){
   radsl<-radiuslow[a]*matrix(1,n,1)  
   uni<-rec2rec(levels[a,],radsl,dendat)
   volslow[a,]<-ccenuni(uni)/maslow[a]
}
vols<-(volshig+volslow)/2
return(vols)
#return(volslow,volshig)
}








