rec2rec<-function(levset,radius,dendat){
#Changes the representation ov levset
#
#"levset" is a n-vector of 0:s an 1:s, indicates the rows of dendat
#radius is n-vector
#dendat is n*d-matrix
#
#Returns setnum*(2*d)-matrix, where setnum is number of sets
#in the union forming the levelset (number of 1:s in levset)
#
len<-length(levset)          #len=n
setnum<-sum(levset)          #number of sets forming the union
levset2<-matrix(0,setnum,1) 
#change the representation of levset from 1/0 to row-numbers of dendat
j<-1
for (i in 1:len){
    if (levset[i]==1){
       levset2[j]<-i
       j<-j+1
    }
}
#Find the endpoints of rectangles
d<-length(dendat[1,])        #d is the dimension of the space
tulos<-matrix(0,setnum,2*d)
for (i in 1:setnum){
   ind<-levset2[i]        #index to row of dendat
   keski<-dendat[ind,]
   r<-radius[ind]
   for (j in 1:d){  
      tulos[i,2*j-1]<-keski[j]-r
      tulos[i,2*j]<-keski[j]+r
   }
}
return(tulos)
}
