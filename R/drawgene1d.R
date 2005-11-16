drawgene1d<-function(values,recs){
#piirtaa 1-ulotteisessa tapauksessa pal. vakion funktion
#
#ohje:   
#bi<-kerbin(dendat,h,binlkm,Q,cfre=FALSE)
#values<-bi$values
#values<-quanti(values,Q,exp(1))
#recs<-bi$recs    
#
recnum<-dim(recs)[1]
x<-matrix(0,2*recnum,1)
for (i in 1:recnum){
   x[i]<-recs[i,1]
   x[recnum+i]<-recs[i,2]
}
y<-matrix(0,2*recnum,1)
y[1:recnum]<-values
y[(recnum+1):(2*recnum)]<-values
plot(x,y,type="h",xlab="",ylab="")
#
for (i in 1:recnum){
    segments(x[i],y[i],x[recnum+i],y[recnum+i])
}
}
