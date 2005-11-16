profgra1<-function(dendat,h,levnum,radius,cfre=FALSE,cvol=FALSE,ccen=FALSE){
#
#dendat<-matrix(rnorm(20),10)
#
eval<-kerest(dendat,h,dendat)       #lask ydinestim arvot hav.pist
#
estim<-makeker1(eval,levnum,radius)  #ydinestim pal vak muodossa
lsets<-estim$lsets
levels<-estim$levels
radius<-estim$radius
#
library(mva)
dismat<-dist(dendat,method="euclidean")
#
alkublokki<-100
blokki<-20
links<-toucbal1(dismat,lsets,radius,alkublokki,blokki)
#
alkublokki2<-200
blokki2<-50
dentree<-decomgra1(lsets,levels,links,radius,alkublokki2,blokki2)
seplsets<-dentree$lsets
sepval<-dentree$levels
parents<-dentree$parents
newradius<-dentree$radius
#
if (cfre) nodefrek<-cfrekv2(seplsets) else nodefrek<-NULL
#
if (cvol){
  vols<-cvolgra1(seplsets,newradius,dendat)
  maslow<-vols$vollow
  mashig<-vols$volhig
  volum<-t((maslow+mashig)/2)   
  kerroin<-cinte(sepval,volum,parents)
  sepval<-sepval/kerroin
}
else volum<-NULL
#
if (ccen && cvol) centers<-ccengra1(seplsets,newradius,dendat,maslow,mashig) 
else centers<-NULL
#                                                                     
return(list(parent=parents,level=sepval,nodefrek=nodefrek,
volume=volum,center=centers))
#values: normeeratut arvot
#nodefrek: kunkin solmun frekvenssi
}       










