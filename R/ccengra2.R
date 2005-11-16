ccengra2<-function(seplsets,newradius,dendat,lowlow,lowhig,higlow,highig,approx=T){
#Laskee joukolle tasojoukon osia painopisteen. 
#
#"seplsets" is tasolkm*n-matriisi, 0:a ja 1:a
#"newradius" is n*tasolkm-matrix
#"dendat" is n*d-matrix
#maslow, mashig are tasolkm-vectors: 
#masses of sets for radius=radiushig, radius=radiuslow correspondingly
#
#Returns a tasolkm*d-matrix: center of mass for each set.
#
#calls: ccenuni
#
n<-length(dendat[,1]) 
radiushig<-newradius
radiuslow<-newradius/sqrt(2)
tasolkm<-length(seplsets[,1])     #seplsets:n rivien maara
d<-length(dendat[1,])
vollowlow<-matrix(0,tasolkm,d)
vollowhig<-matrix(0,tasolkm,d)
volhiglow<-matrix(0,tasolkm,d)
volhighig<-matrix(0,tasolkm,d)
for (a in 1:tasolkm){
  radsh<-radiushig[,a]
  uni<-rec2rec(seplsets[a,],radsh,dendat)
  if (approx){
     appu<-ccenuniapp(uni)
     volhiglow[a,]<-appu$low/higlow[a]
     volhighig[a,]<-appu$hig/highig[a]
  }
  else vollowlow[a]<-ccenuni(runi)/mashig[a]
}
for (a in 1:tasolkm){
   radsl<-radiuslow[,a]
   uni<-rec2rec(seplsets[a,],radsl,dendat)
   if (approx){
     appu<-ccenuniapp(uni)
     vollowlow[a,]<-appu$low/lowlow[a]
     vollowhig[a,]<-appu$hig/lowhig[a]
    }
    else vollowlow[a,]<-ccenuni(uni)/lowlow[a]
}
return(list(vollowlow=vollowlow,vollowhig=vollowhig,volhiglow=volhiglow,
volhighig=volhighig))
}








