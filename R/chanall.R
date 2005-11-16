chanall<-function(lsets){
#
if (dim(t(lsets))[1]==1) levnum<-1 else levnum<-length(lsets[,1]) 
len<-length(lsets[1,])
#
ma<-0
lsetsco<-matrix(NA,levnum,len)
maara<-matrix(0,levnum,1)   #need for drawing gra
for (i in 1:levnum){
    maara[i]<-sum(lsets[i,]) 
    j<-1
    for (k in 1:len){
        if (lsets[i,k]==1){
           lsetsco[i,j]<-k
           j<-j+1
        }
    ma<-max((j-1),ma)
    }
}
#lsetsco<-lsetsco[1:ma,]
return(lsetsco=lsetsco,maara=maara)
}

