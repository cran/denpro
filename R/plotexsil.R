plotexsil<-function(ltree,n){
#
level<-1
crit<-rep(0,d)
#
vecplu<-prof2vecs(ltree,level,n,crit)
vecs<-vecplu$vecs
depths<-vecplu$depths
motes<-vecplu$motes
#
return(list(vecs=vecs,depths=t(depths),motes=t(motes)))
}
     
