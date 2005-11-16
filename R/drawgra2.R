drawgra2<-function(dendat,h,levnum,ks,plkm){
#
#ohje:   dendat<-matrix(rnorm(20),10)
#        dr<-drawgra2(dendat,h=1,levnum=3,ks=c(3,2,1),plkm=30)
#        persp(dr$x,dr$y,dr$z,theta=130,phi=30)  
#
eval<-kerest(dendat,h,dendat)       #lask ydinestim arvot hav.pist
#
df<-dfind(dendat,ks)
close<-df$close  
radius<-df$radius  
#
estim<-makeker2(eval,levnum,ks,close) 
lsets<-estim$lsets
levels<-estim$levels          
remai<-length(lsets[,1])   #top lev.sets may be cutted
radius<-radius[,1:remai]  
M<-estim$max
#
muunnos<-chanall(lsets)
lpoints2<-muunnos$lsetsco
amount<-muunnos$maara
#
mrad<-max(radius)                   #maximal radius
xmax<-max(dendat[,1])+mrad
xmin<-min(dendat[,1])-mrad
ymax<-max(dendat[,2])+mrad
ymin<-min(dendat[,2])-mrad
xaskel<-(xmax-xmin)/(plkm-1)
yaskel<-(ymax-ymin)/(plkm-1)
x<-seq(xmin,xmax,xaskel) 
y<-seq(ymin,ymax,yaskel)
z<-matrix(0,length(x),length(y))
#
lnum<-length(lsets[,1])   #number of lsets
for (l in 1:lnum){
   for (m in 1:amount[l]){
      obsind<-lpoints2[l,m]
      xalku<-dendat[obsind,1]-radius[obsind,l]
      yalku<-dendat[obsind,2]-radius[obsind,l]
      xloppu<-dendat[obsind,1]+radius[obsind,l]
      yloppu<-dendat[obsind,2]+radius[obsind,l]
      xkeski<-(xalku+xloppu)/2
      ykeski<-(yalku+yloppu)/2 
      xalkui<-floor((xalku-xmin)/xaskel)+1
      xloppui<-floor((xloppu-xmin)/xaskel)+1
      yalkui<-floor((yalku-ymin)/yaskel)+1
      yloppui<-floor((yloppu-ymin)/yaskel)+1
      for (i in xalkui:xloppui){
        for (j in yalkui:yloppui){
          if ((x[i]-xkeski)^2+(y[j]-ykeski)^2<=radius[l]^2) z[i,j]<-levels[l]
        }
      }
   }
}
return(list(x=x,y=y,z=z))
#return(list(x=x,y=y,z=z,tieto=levels,max=M)) 
}











