scauniv<-function(dendat,hvec,binvec,Q){
#
hnum<-length(hvec)
for (i in 1:hnum){   
    h<-hvec[i]
    binlkm<-binvec[i]
    plkm<-50
    dr<-drawdisk1d(dendat,h,binlkm,Q,plkm)
    plot(dr$x,dr$z,type="b")
}
}
