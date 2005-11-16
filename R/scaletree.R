scaletree<-function(bt,suppo,parvec=NULL,makeplot=F,treeseq=NULL){
#
d<-length(bt$step)
#d<-dim(dendat)[2]
#n<-dim(dendat)[1]
#
if (is.null(treeseq)){
  treeseq<-prune(bt)
}
#suppo<-supp(dendat)
if (is.null(parvec)){
   parvec<-treeseq$leafs
}
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
     leafnum<-parvec[i]
     tree<-picktree(treeseq,leafnum)
     #
     pv<-partition(tree,suppo)
     recs<-pv$recs
     values<-pv$values  
     #
     profile<-profgene(values,recs,frekv=F,cvol=T,ccen=T,cfre=F)
     #
     level<-leafnum
     vecplu<-prof2vecs(profile,level,n,crit)
     vecs<-vecplu$vecs
     depths<-vecplu$depths
     motes<-vecplu$motes
     mlabel<-vecplu$mlabel
     vecnum<-length(depths)
     smoot<-matrix(level,vecnum,1)
     #
     # concatenate to big's
     #
     veclkmold<-veclkm
     veclkm<-veclkm+vecnum
     if (veclkmold==0){   
        bigvecs<-vecs
        bigdepths<-depths
        bigmotes<-motes
        bigmlabel<-mlabel
        bigsmoot<-smoot
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
       templab<-matrix(0,veclkm,1)
       templab[1:veclkmold]<-bigmlabel
       templab[(veclkmold+1):veclkm]<-mlabel
       bigmlabel<-templab
       #
       tempsmoo<-matrix(0,veclkm,1)
       tempsmoo[1:veclkmold]<-bigsmoot
       tempsmoo[(veclkmold+1):veclkm]<-smoot
       bigsmoot<-tempsmoo
     }
}    
if (makeplot==T) plotvecs(bigvecs,segme=T,bigdepths) 
#
return(list(bigvecs=bigvecs,bigdepths=t(bigdepths),motes=t(bigmotes),mlabel=t(bigmlabel),smoot=t(bigsmoot)))
}







