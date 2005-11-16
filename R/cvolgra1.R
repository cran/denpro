cvolgra1<-function(levels,radius,dendat){
#Laskee joukolle tasojoukon osia voluumin.
#
#"levels" is tasolkm*n-matriisi, 0:a ja 1:a
#"radius" is tasolkm-vector
#"dendat" is n*d-matrix
#
#Returns a tasolkm-vector: volume for each set.
#            
n<-length(dendat[,1])
radiushig<-radius
radiuslow<-radius/sqrt(2)
tasolkm<-length(levels[,1])     #levels:n rivien maara
volshig<-matrix(0,tasolkm,1)
volslow<-matrix(0,tasolkm,1)
for (a in 1:tasolkm){
  radsh<-radiushig[a]*matrix(1,n,1)
  uni<-rec2rec(levels[a,],radsh,dendat)
  volshig[a]<-til(uni)
}
for (a in 1:tasolkm){
   radsl<-radiuslow[a]*matrix(1,n,1)
   uni<-rec2rec(levels[a,],radsl,dendat)
   volslow[a]<-til(uni)
}
#vols<-(volshig+volslow)/2
#return(t(vols))
return(list(vollow=volslow,volhig=volshig))
}             
