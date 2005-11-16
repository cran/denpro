profdisk<-function(dendat,h,binlkm,Q,cvol=TRUE,ccen=TRUE,cfre=FALSE){
#esim. dendat<-matrix(rnorm(20),10) on 10*2 matriisi
#
bi<-kerbin(dendat,h,binlkm,Q,cfre)
values<-bi$values
recs<-bi$recs
frekv<-bi$frek 
#
pr<-profgene(values=values,recs=recs,frekv=frekv,cvol=cvol,ccen=ccen,
cfre=cfre)
#
return(list(parents=pr$parents,levels=pr$levels,invalues=pr$invalues,
volumes=pr$volumes,centers=pr$centers,nodefrek=pr$nodefrek))
}









