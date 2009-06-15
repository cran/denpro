paracoor<-function(dendat,xmargin=0.1,
paletti=matrix("black",dim(dendat)[1],1),noadd=TRUE,verti=NULL,cex.axis=1)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
ylim<-c(min(dendat),max(dendat))

if (noadd)
plot(x="",y="",
xlim=c(1-xmargin,d+xmargin),ylim=ylim,
xlab="",ylab="",xaxt="n",cex.axis=cex.axis)

for (i in 1:n){
      points(dendat[i,],col=paletti[i])
      for (j in 1:(d-1)) segments(j,dendat[i,j],j+1,dendat[i,j+1],
      col=paletti[i])
}

if (!is.null(verti)) segments(verti,ylim[1],verti,ylim[2])

}

