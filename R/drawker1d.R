drawker1d<-function(dendat,h,plkm,kernel="normal"){
#
sade<-2
xmax<-max(dendat)+sade
xmin<-min(dendat)-sade
xaskel<-(xmax-xmin)/(plkm-1)
x<-seq(xmin,xmax,xaskel) 
if (kernel=="normal"){
    z<-kerest1d(x,h,dendat) 
}
#
return(list(x=x,z=z))
}

