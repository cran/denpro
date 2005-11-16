scaledisk<-function(dendat,parvec,binvec,Q,makeplot=F){
#
#hved is hnum-vector, h:s from small to large
#
d<-dim(dendat)[2]
n<-dim(dendat)[1]
#
parnum<-length(parvec)
#
veclkm<-0
if (d==1){
  crit<-min(dendat)
}
else{
  crit<-rep(0,d)
}
for (i in 1:parnum){ 
     h<-parvec[i]
     binlkm<-binvec[i]
     profile<-profdisk(dendat,h,binlkm,Q,cvol=TRUE,ccen=TRUE,cfre=FALSE)
     #
     level<-log(h,base=10)
     vecplu<-prof2vecs(profile,level,n,crit)
     vecs<-vecplu$vecs
     depths<-vecplu$depths
     motes<-vecplu$motes
     vecnum<-length(depths)
     #
     # concatenate to big's
     #
     veclkmold<-veclkm
     veclkm<-veclkm+vecnum
     if (veclkmold==0){   
        bigvecs<-vecs
        bigdepths<-depths
        bigmotes<-motes
     }
     else{
       temp<-matrix(0,veclkm,4)
       temp[1:veclkmold,]<-bigvecs
       temp[(veclkmold+1):veclkm,]<-vecs
       bigvecs<-temp
       #
       tempdep<-matrix(0,veclkm,1)
       tempdep[1:veclkmold]<-bigdepths
       tempdep[(veclkmold+1):veclkm]<-depths
       bigdepths<-tempdep
       #
       tempmoo<-matrix(0,veclkm,1)
       tempmoo[1:veclkmold]<-bigmotes
       tempmoo[(veclkmold+1):veclkm]<-motes
       bigmotes<-tempmoo
       #
     }
}    
if (makeplot==T) plotvecs(bigvecs,segme=T,bigdepths) 
#
return(list(bigvecs=bigvecs,bigdepths=t(bigdepths),motes=t(bigmotes)))
}







