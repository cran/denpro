profgran<-function(dendat,h,levnum,ks,cfre=FALSE,cvol=TRUE,ccen=TRUE){
#
#dendat<-matrix(rnorm(20),10)
#
eval<-kerest(dendat,h,dendat)       #lask ydinestim arvot hav.pist
#
df<-dfind(dendat,ks)
close<-df$close
radius<-df$radius
dismat<-df$dismat      #library(mva)#dismat<-dist(dendat,method="euclidean")
#
estim<-makeker2(eval,levnum,ks,close)  #ydinestim pal vak muodossa
lsets<-estim$lsets
levels<-estim$levels
remai<-length(lsets[,1])    #top lev.sets may be cutted
radius<-radius[,1:remai]
#
alkublokki<-100
blokki<-20
links<-toucbal2(dismat,lsets,radius,alkublokki,blokki)
#
alkublokki2<-200
blokki2<-50
dentree<-decomgra2(lsets,levels,links,radius,alkublokki2,blokki2)
seplsets<-dentree$lsets
sepval<-dentree$levels
parents<-dentree$parents
newradius<-dentree$radius
#
if (cfre) nodefrek<-cfrekv2(seplsets) else nodefrek<-NULL
#
if (ccen) cvol<-T
if (cvol){
    vols<-cvolgra2(seplsets,newradius,dendat,approx=T)
    lowlow<-vols$lowlow
    lowhig<-vols$lowhig
    higlow<-vols$higlow
    highig<-vols$highig
    volum<-t((lowlow+lowhig+higlow+highig)/4)   
    kerroin<-cinte(sepval,volum,parents)
    sepval<-sepval/kerroin
}
else volum<-NULL
#
if (ccen){
   mcenters<-ccengra2(seplsets,newradius,dendat,lowlow,lowhig,higlow,highig,
                     approx=T) 
   vollowlow<-mcenters$vollowlow
   vollowhig<-mcenters$vollowhig
   volhiglow<-mcenters$volhiglow
   volhighig<-mcenters$volhighig
   centers<-t((vollowlow+vollowhig+volhiglow+volhighig)/4)   
}
else{
  centers<-NULL
}
#                                                                     
return(list(parent=parents,level=sepval,nodefrek=nodefrek,
volume=volum,center=centers))
#values: normeeratut arvot
#nodefrek: kunkin solmun frekvenssi
}       










