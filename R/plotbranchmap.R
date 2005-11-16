plotbranchmap<-function(bm,phi=55,theta=30)
{

persp(x=bm$level,y=bm$h,z=bm$z, 
xlab="level",ylab="h",zlab="excess mass",
ticktype="detailed",
col=bm$col,phi=phi,theta=theta) 

}


