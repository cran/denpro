paracoor<-function(dendat,xmargin=0.1,
paletti=matrix("black",dim(dendat)[1],1))
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

plot(x="",y="",
xlim=c(1-xmargin,d+xmargin),ylim=c(min(dendat),max(dendat)),
xlab="",ylab="",xaxt="n")

for (i in 1:n){
      points(dendat[i,],col=paletti[i])
      for (j in 1:(d-1)) segments(j,dendat[i,j],j+1,dendat[i,j+1],
      col=paletti[i])
}

}

