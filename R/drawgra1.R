drawgra1<-function(dendat,h,levnum,radius,plkm){
#
#ohje:   dendat<-matrix(rnorm(20),10)
#        dr<-drawgra1(dendat,h=1,levnum=3,radius=c(1,.5,.3),plk=30)
#        persp(dr$x,dr$y,dr$z,theta=130,phi=30)  
#
eval<-kerest(dendat,h,dendat)       #lask ydinestim arvot hav.pist
#
estim<-makeker1(eval,levnum,radius)  #ydinestim pal vak muodossa
lsets<-estim$lsets
levels<-estim$levels
radius<-estim$radius           
M<-estim$max
#
muunnos<-chanall(lsets)
lpoints2<-muunnos$lsetsco
amount<-muunnos$maara
#
xmax<-max(dendat[,1])+radius[1]
xmin<-min(dendat[,1])-radius[1]
ymax<-max(dendat[,2])+radius[1]
ymin<-min(dendat[,2])-radius[1]
xaskel<-(xmax-xmin)/(plkm-1)
yaskel<-(ymax-ymin)/(plkm-1)
x<-seq(xmin,xmax,xaskel) 
y<-seq(ymin,ymax,yaskel)
z<-matrix(0,length(x),length(y))
#
lnum<-length(lsets[,1])   #number of lsets
for (l in 1:lnum){
   for (m in 1:amount[l]){
      xalku<-dendat[lpoints2[l,m],1]-radius[l]
      yalku<-dendat[lpoints2[l,m],2]-radius[l]
      xloppu<-dendat[lpoints2[l,m],1]+radius[l]
      yloppu<-dendat[lpoints2[l,m],2]+radius[l]
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











