drawdisk<-function(dendat,h,binlkm,Q,plkm){
#piirtaa 2-ulotteisessa tapauksessa diskretoidun ydinestimaattorin kuvaajan 
#
#plkm on kuvaajan hilan pisteiden lkm
#
#ohje:   dendat<-matrix(rnorm(20),10) 
#        dr<-drawbink(dendat,h=1,binlkm=5,Q=3,plk=30)
#        persp(dr$x,dr$y,dr$z,theta=130,phi=30)
#
kanta<-support(dendat,epsi=0)
kanta1<-kanta[1,1]-h
kanta2<-kanta[1,2]+h
kanta3<-kanta[2,1]-h
kanta4<-kanta[2,2]+h
pit1<-(kanta2-kanta1)/plkm
pit2<-(kanta4-kanta3)/plkm
x<-kanta1+c(0:plkm)*pit1
y<-kanta3+c(0:plkm)*pit2
#
bi<-kerbin(dendat,h,binlkm,Q,cfre=FALSE)
values<-bi$values
values<-quanti(values,Q,exp(1))
recs<-bi$recs    
#
ep2<-.01
xdim<-length(x)
ydim<-length(y)
arvot<-matrix(0,xdim,ydim)
lnum<-length(recs[,1])
#
i<-1
while (i<=xdim){
  j<-1
  while (j<=ydim){
    k<-1
    while (k<=lnum){
      if ((x[i]<=recs[k,2]+ep2) && (x[i]>=recs[k,1]-ep2) &&
          (y[j]<=recs[k,4]+ep2) && (y[j]>=recs[k,3]-ep2))
         arvot[i,j]<-values[k]
    k<-k+1
    }
  j<-j+1
  }
i<-i+1
}
return(list(x=x,y=y,z=arvot))
#persp(x,y,arvot)
}







