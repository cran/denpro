paracoor<-function(dendat)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

plot(x="",y="",
xlim=c(0.5,d+0.5),ylim=c(min(dendat),max(dendat)),
xlab="",ylab="",xaxt="n")

for (i in 1:n){
      points(dendat[i,]) #,col="red")
      for (j in 1:(d-1)) segments(j,dendat[i,j],j+1,dendat[i,j+1])
}

}

