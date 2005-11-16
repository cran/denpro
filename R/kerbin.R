kerbin<-function(dendat,h,binlkm,Q,cfre=FALSE){
#Constructs a binned kernel estimator
#
#if (dim(t(dendat))[1]==1) d<-1 else d<-length(dendat[1,]) 
d<-dim(dendat)[2]
if (d==1){
   palvak<-makedk1d(dendat,h,binlkm)
}
else{
   palvak<-makedk(dendat,h,binlkm)
}
values<-palvak$values
recs<-palvak$recs
#
values<-quanti(values,Q,exp(1)) #quantisising the kernel estimate
pvak<-rmzeros(values,recs)      #removing zeros
values<-pvak$values
recs<-pvak$recs
#
return(values=values,recs=recs,frek=NULL)
}   
