rmzeros<-function(values,recs){
#
itemnum<-length(values)
d2<-length(recs[1,])
newvalues<-matrix(0,itemnum,1)
newrecs<-matrix(0,itemnum,d2)
ind<-1
for (i in 1:itemnum){
  if (values[i]>0){
     newvalues[ind]<-values[i]
     newrecs[ind,]<-recs[i,]
     ind<-ind+1
  }
}
newvalues<-newvalues[1:(ind-1)]
newrecs<-newrecs[1:(ind-1),]
return(list(values=newvalues,recs=newrecs))
}






