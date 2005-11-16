plotgrapre<-function(dendat,h,levnum,ks){
#
#ALKU OF profgra2
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
newradius<-dentree$radius
#                                 
# END OF profgra2 SEGMENT
#
return(list(seplsets=seplsets,newradius=newradius,sepval=sepval))
}
